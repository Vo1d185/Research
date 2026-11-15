"""Isoform inventory of the LR panel genes (UniProt, mouse). Counts and flags only; no analysis.

For every gene of a prepared stage-2 input (ligands and receptors; complexes split into subunits):

  * reviewed mouse UniProt entries of the gene: entries whose primary gene name is the gene (proteins.search),
    plus every reviewed mouse entry cross-referenced to the same MGI gene as one of those entries
    (query xref:"mgi-MGI:..."). Several reviewed entries of one gene (e.g. Nrxn1/2/3 alpha and beta, separate
    entries from alternative promoters) are treated as different isoforms of that gene;
  * isoforms: ALTERNATIVE PRODUCTS (canonical + alternative products) of each entry; for the gene, the union of
    isoform identifiers over its entries (an entry without the annotation counts as one isoform, <accession>-1);
  * variable regions, two sources:
      - alternative_sequence: 'Alternative sequence' features of an entry (that entry's canonical coordinates);
      - between_entries: for a gene with several entries, each entry is compared with the gene's reference entry
        (the longest; ties by accession). Identical co-linear stretches of >= MIN_BLOCK residues
        (difflib.SequenceMatcher, autojunk off) are shared; every reference stretch between shared blocks is a
        region (absent from or replaced in the other entry), and residues of the other entry with no reference
        counterpart between two shared blocks are an insertion at that reference position.

Each variable region is flagged (in the coordinates of the entry it is described on) when it overlaps
  * the extracellular region: 'Topological domain' features described as Extracellular; for an entry
    without topological domains whose subcellular location includes Secreted, the chain after the
    signal peptide (the whole chain if no signal peptide is annotated);
  * an annotated domain ('Domain' features);
  * an annotated binding region ('Region' features whose description contains 'bind' or 'interact')
    or a 'Binding site' feature.

A gene is a candidate when (a) it has >= 2 isoforms in one entry and an alternative_sequence region with any flag
(inventory 01 rule), or (b) it has >= 2 reviewed entries and a between_entries region with any flag.
Nrxn1 alpha/beta is the first row of candidates.csv as the known example.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import difflib
import hashlib
import json
from pathlib import Path
import time
from urllib.error import HTTPError
from urllib.parse import quote

import pandas as pd

import proteins
from research import check, launch_info

KNOWN = dict(gene='Nrxn1', note='known example: Nrxn1 alpha (Q9CS84) and beta (P0DI97) are separate reviewed '
             'entries (alternative promoters); see HANDOFF.md (beta used in stage 1)')
MIN_BLOCK = 10          # identical co-linear stretch (residues) counted as shared between two entries of a gene
FLAGS = ['overlaps_extracellular', 'overlaps_domain', 'overlaps_binding']


def polite(fetch, tries=8):
    """Retry on HTTP 429/503 (UniProt rate limit), honouring Retry-After, with exponential backoff."""
    def call(url):
        for attempt in range(tries):
            try:
                return fetch(url)
            except HTTPError as error:
                if error.code not in (429, 503) or attempt == tries - 1:
                    raise
                wait = error.headers.get('Retry-After') if error.headers else None
                time.sleep(float(wait) if wait and str(wait).isdigit() else 2.0 * 2 ** attempt)
    return call


def panel_genes(input_dir):
    """gene -> (roles, LR pairs) from lr.csv and entities.csv of a prepared stage-2 input."""
    input_dir = Path(input_dir)
    lr = pd.read_csv(input_dir/'lr.csv', dtype=str)
    entities = pd.read_csv(input_dir/'entities.csv', dtype=str)
    members = {r.entity: r.subunits.split(';') for r in entities.itertuples(index=False)}
    genes = {}
    for pair in lr.itertuples(index=False):
        for role, entity in [('ligand', pair.ligand), ('receptor', pair.receptor)]:
            for gene in members.get(entity, [entity]):
                g = genes.setdefault(gene, dict(roles=set(), pairs=set(), types=set()))
                g['roles'].add(role)
                g['pairs'].add(f'{pair.ligand}->{pair.receptor}')
                g['types'].add(pair.type)
    return genes, lr


def mgi_ids(entry):
    return sorted({r.get('id') for r in entry.get('uniProtKBCrossReferences', [])
                   if r.get('database') == 'MGI' and r.get('id')})


def xref_search(mgi_id, fetch):
    """Reviewed mouse entries cross-referenced to one MGI gene."""
    query = quote(f'(xref:"mgi-{mgi_id}") AND (organism_id:{proteins.MOUSE_TAXON}) AND (reviewed:true)')
    found, headers = fetch(proteins.SEARCH.format(query=query))
    hits = []
    for entry in found.get('results', []):
        fields = proteins.entry_fields(entry)
        if fields['taxon_id'] == proteins.MOUSE_TAXON and fields['reviewed'] and mgi_id in mgi_ids(entry):
            hits.append((entry, fields))
    return hits, headers


def gene_entries(gene, fetch):
    """[(entry, fields, found_by)], MGI ids, UniProt releases seen."""
    hits, headers = proteins.search(gene, fetch)
    found = {f['accession']: [e, f, {'primary_name'}] for e, f in hits}
    ids = sorted({m for e, _ in hits for m in mgi_ids(e)})
    releases = {headers.get('x-uniprot-release')}
    for mgi in ids:
        more, h = xref_search(mgi, fetch)
        releases.add(h.get('x-uniprot-release'))
        for e, f in more:
            found.setdefault(f['accession'], [e, f, set()])[2].add('mgi_xref')
    return [(e, f, ';'.join(sorted(s))) for e, f, s in found.values()], ids, releases


def overlaps(a, ranges):
    return [r for r in ranges if a[0] <= r[1] and r[0] <= a[1]]


def labelled(entry, kind, test=None):
    """(start, end, description) of features of one type (1-based, canonical coordinates)."""
    out = []
    for f in entry.get('features', []):
        if f.get('type') != kind:
            continue
        description = str(f.get('description', ''))
        if test and not test(description):
            continue
        start = f.get('location', {}).get('start', {}).get('value')
        end = f.get('location', {}).get('end', {}).get('value')
        if isinstance(start, int) and isinstance(end, int) and 1 <= start <= end:
            out.append((start, end, description))
    return out


def extracellular(entry, length):
    topo = labelled(entry, 'Topological domain')
    if topo:
        return [(s, e) for s, e, d in topo if 'extracellular' in d.lower()], 'topological_domain'
    locations = [loc.get('location', {}).get('value', '') for c in entry.get('comments', [])
                 if c.get('commentType') == 'SUBCELLULAR LOCATION' for loc in c.get('subcellularLocations', [])]
    if any('secreted' in str(v).lower() for v in locations):
        signal = labelled(entry, 'Signal')
        start = max(e for _, e, _ in signal) + 1 if signal else 1
        return [(start, length)], 'secreted_after_signal' if signal else 'secreted_whole_chain'
    return [], 'none_annotated'


def isoforms(entry):
    for c in entry.get('comments', []):
        if c.get('commentType') == 'ALTERNATIVE PRODUCTS':
            return [dict(name=i.get('name', {}).get('value'), ids=';'.join(i.get('isoformIds', [])),
                         status=i.get('isoformSequenceStatus')) for i in c.get('isoforms', [])]
    return []


def isoform_ids(entry, accession):
    ids = {i for iso in isoforms(entry) for i in iso['ids'].split(';') if i}
    return ids or {f'{accession}-1'}


def annotations(entry, length):
    ec, basis = extracellular(entry, length)
    domains = labelled(entry, 'Domain')
    binding = labelled(entry, 'Region', lambda d: 'bind' in d.lower() or 'interact' in d.lower()) + \
        labelled(entry, 'Binding site')
    return ec, basis, domains, binding


def flag(region, ec, domains, binding):
    return dict(overlaps_extracellular=bool(overlaps(region, ec)),
                overlaps_domain=bool(overlaps(region, [(s, e) for s, e, _ in domains])),
                overlaps_binding=bool(overlaps(region, [(s, e) for s, e, _ in binding])),
                domains=';'.join(sorted({d for s, e, d in domains if overlaps(region, [(s, e)])})),
                binding=';'.join(sorted({d for s, e, d in binding if overlaps(region, [(s, e)])})))


def entry_differences(reference, other, min_block=MIN_BLOCK):
    """Regions of `reference` not shared with `other` (1-based, inclusive reference coordinates), and the number
    of shared residues. Shared = identical co-linear blocks of >= min_block residues. An insertion in `other`
    between reference residues k and k+1 is reported as the region (k, k+1) (k = 0 or len: the end residue)."""
    blocks = [b for b in difflib.SequenceMatcher(None, reference, other, autojunk=False).get_matching_blocks()
              if b.size >= min_block]
    regions, ra, oa = [], 0, 0
    for b in blocks + [difflib.Match(len(reference), len(other), 0)]:
        ref_gap, other_gap = b.a - ra, b.b - oa
        if ref_gap > 0 or other_gap > 0:
            if ref_gap > 0:
                start, end = ra + 1, b.a
            else:
                start, end = max(ra, 1), min(ra + 1, len(reference))
            kind = 'absent' if other_gap == 0 else 'inserted' if ref_gap == 0 else 'replaced'
            regions.append(dict(start=start, end=end, kind=kind, reference_residues=ref_gap,
                                other_start=oa + 1 if other_gap else None, other_end=b.b if other_gap else None,
                                other_residues=other_gap))
        ra, oa = b.a + b.size, b.b + b.size
    return regions, sum(b.size for b in blocks)


def describe(gene, entry, fields, found_by):
    length = len(fields['sequence'])
    ec, basis, domains, binding = annotations(entry, length)
    variants = []
    for start, end, description in labelled(entry, 'Alternative sequence'):
        variants.append(dict(gene=gene, source='alternative_sequence', accession=fields['accession'],
                             compared_accession=None, start=start, end=end, description=description,
                             **flag((start, end), ec, domains, binding)))
    iso = isoforms(entry)
    return dict(gene=gene, accession=fields['accession'], uniprot_id=fields['uniprot_id'],
                protein_name=fields['protein_name'], gene_primary=fields['gene_primary'], found_by=found_by,
                mgi_ids=';'.join(mgi_ids(entry)), length=length, isoforms=max(len(iso), 1),
                isoform_names=';'.join(str(i['name']) for i in iso), variable_regions=len(variants),
                extracellular_basis=basis, extracellular_ranges=';'.join(f'{s}-{e}' for s, e in ec),
                entry_version=fields['entry_version'], sequence_sha256=proteins.sha256_text(fields['sequence'])), variants


def between_entries(gene, hits):
    """Variable regions between the reference entry (longest; ties by accession) and each other entry."""
    if len(hits) < 2:
        return [], None
    ref_entry, ref = sorted(((e, f) for e, f, _ in hits), key=lambda h: (-len(h[1]['sequence']), h[1]['accession']))[0]
    ec, _, domains, binding = annotations(ref_entry, len(ref['sequence']))
    out = []
    for _, f, _ in hits:
        if f['accession'] == ref['accession']:
            continue
        regions, shared = entry_differences(ref['sequence'], f['sequence'])
        for r in regions:
            where = f"{f['accession']} {r['other_start']}-{r['other_end']}" if r['other_residues'] else f"not in {f['accession']}"
            out.append(dict(gene=gene, source='between_entries', accession=ref['accession'],
                            compared_accession=f['accession'], start=r['start'], end=r['end'],
                            description=f"{r['kind']}: {ref['accession']} {r['start']}-{r['end']} "
                                        f"({r['reference_residues']} aa) -> {where} ({r['other_residues']} aa); "
                                        f"shared {shared} aa in blocks >= {MIN_BLOCK}",
                            kind=r['kind'], reference_residues=r['reference_residues'],
                            other_residues=r['other_residues'], shared_residues=shared,
                            **flag((r['start'], r['end']), ec, domains, binding)))
    return out, ref['accession']


def any_flag(variants, gene, source):
    v = variants[(variants.gene == gene) & (variants.source == source)] if len(variants) else variants
    return {c: bool(v[c].any()) if len(v) else False for c in FLAGS}


def run(args, fetch=None):
    launch = launch_info()
    fetch = fetch or polite(proteins.http_json)
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    genes, lr = panel_genes(args.input)
    names = sorted(genes)
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(lambda g: (g,) + gene_entries(g, fetch), names))
    entries, variants, raw, releases, per_gene = [], [], {}, set(), []
    for gene, hits, ids, rel in results:
        releases |= rel
        if not hits:
            entries.append(dict(gene=gene, accession=None, protein_name='no reviewed mouse entry with this primary gene name'))
        merged = set()
        for entry, fields, found_by in hits:
            row, var = describe(gene, entry, fields, found_by)
            entries.append(row)
            variants.extend(var)
            merged |= isoform_ids(entry, fields['accession'])
            raw[fields['accession']] = proteins.sha256_text(json.dumps(entry, sort_keys=True))
        var, reference = between_entries(gene, hits)
        variants.extend(var)
        per_gene.append(dict(gene=gene, reviewed_entries=len(hits),
                             accessions=';'.join(f['accession'] for _, f, _ in hits),
                             entries_by_mgi_only=sum(fb == 'mgi_xref' for _, _, fb in hits),
                             mgi_ids=';'.join(ids), reference_entry=reference, merged_isoforms=len(merged),
                             max_isoforms=max([r['isoforms'] for r in entries if r['gene'] == gene and r.get('accession')],
                                              default=None)))
    entries, variants, per_gene = pd.DataFrame(entries), pd.DataFrame(variants), pd.DataFrame(per_gene)
    if variants.empty:
        variants = pd.DataFrame(columns=['gene', 'source', 'accession', 'compared_accession', 'start', 'end'] + FLAGS)
    alt = variants[variants.source == 'alternative_sequence']
    flags = alt.groupby(['gene', 'accession'])[FLAGS].any().reset_index()
    entries = entries.merge(flags, on=['gene', 'accession'], how='left')
    for c in FLAGS:
        entries[c] = entries[c].astype('boolean').fillna(False).astype(bool)
    per_gene['variable_regions'] = per_gene.gene.map(lambda g: int((alt.gene == g).sum()))
    per_gene['between_entry_regions'] = per_gene.gene.map(lambda g: int(((variants.gene == g) &
                                                                         (variants.source == 'between_entries')).sum()))
    for source, prefix in [('alternative_sequence', 'any'), ('between_entries', 'between')]:
        f = per_gene.gene.map(lambda g: any_flag(variants, g, source))
        for c, name in zip(FLAGS, ['extracellular', 'domain', 'binding']):
            per_gene[f'{prefix}_{name}'] = f.map(lambda d: d[c])
    rule_a = (per_gene.max_isoforms.fillna(0) >= 2) & (per_gene.any_extracellular | per_gene.any_domain |
                                                       per_gene.any_binding)
    rule_b = (per_gene.reviewed_entries >= 2) & (per_gene.between_extracellular | per_gene.between_domain |
                                                 per_gene.between_binding)
    per_gene['candidate_alternative_sequence'] = rule_a
    per_gene['candidate_several_entries'] = rule_b
    per_gene['candidate'] = rule_a | rule_b
    per_gene['candidate_rule'] = [';'.join(n for n, v in [('alternative_sequence', a), ('several_entries', b)] if v)
                                  for a, b in zip(rule_a, rule_b)]
    per_gene['roles'] = per_gene.gene.map(lambda g: ';'.join(sorted(genes[g]['roles'])))
    per_gene['types'] = per_gene.gene.map(lambda g: ';'.join(sorted(genes[g]['types'])))
    per_gene['n_pairs'] = per_gene.gene.map(lambda g: len(genes[g]['pairs']))
    per_gene['lr_pairs'] = per_gene.gene.map(lambda g: ';'.join(sorted(genes[g]['pairs'])))
    candidates = per_gene[per_gene.candidate].copy()
    known = per_gene[per_gene.gene == KNOWN['gene']].copy()
    if len(known):
        known['note'] = KNOWN['note']
        candidates = pd.concat([known, candidates[candidates.gene != KNOWN['gene']]], ignore_index=True)
    candidates['note'] = candidates.get('note', pd.Series(dtype=str)).fillna('')
    changes = None
    if args.previous:
        prev = pd.read_csv(Path(args.previous)/'genes.csv').set_index('gene')
        changes = per_gene.set_index('gene')[['accessions', 'candidate', 'candidate_rule']].join(
            prev[['accessions', 'candidate']].rename(columns=lambda c: c + '_previous'), how='outer')
        changes = changes[(changes.candidate != changes.candidate_previous) |
                          (changes.accessions.fillna('') != changes.accessions_previous.fillna(''))].reset_index()
    out.mkdir(parents=True)
    per_gene.to_csv(out/'genes.csv', index=False)
    entries.to_csv(out/'entries.csv', index=False)
    variants.to_csv(out/'variable_regions.csv', index=False)
    candidates.to_csv(out/'candidates.csv', index=False)
    if changes is not None:
        changes.to_csv(out/'changes_vs_previous.csv', index=False)
    prep = json.loads((Path(args.input)/'preparation.json').read_text())
    (out/'inventory.json').write_text(json.dumps(dict(status='complete', launch=launch,
        created=datetime.now(timezone.utc).isoformat(timespec='seconds'), input=str(args.input),
        sample=prep.get('sample'), input_sha256={f: hashlib.sha256((Path(args.input)/f).read_bytes()).hexdigest()
                                                 for f in ['lr.csv', 'entities.csv']},
        uniprot_releases=sorted(r for r in releases if r),
        query='proteins.search: (gene_exact:G) AND (organism_id:10090) AND (reviewed:true), primary gene name must '
              'equal G; then for each MGI id of those entries: (xref:"mgi-MGI:...") AND (organism_id:10090) AND '
              '(reviewed:true), entry must carry that MGI id',
        min_block=MIN_BLOCK, previous=args.previous,
        genes=len(per_gene), genes_with_entry=int((per_gene.reviewed_entries > 0).sum()),
        genes_with_several_entries=sorted(per_gene[per_gene.reviewed_entries > 1].gene),
        entries_found_by_mgi_only=sorted(entries[entries.get('found_by', pd.Series(dtype=str)) == 'mgi_xref']
                                         .accession.dropna()),
        candidates=int(per_gene.candidate.sum()),
        candidates_alternative_sequence=int(rule_a.sum()), candidates_several_entries=int(rule_b.sum()),
        candidates_several_entries_only=sorted(per_gene[rule_b & ~rule_a].gene),
        entry_json_sha256=raw), indent=2, default=str))
    print(f'{len(per_gene)} genes, {int(per_gene.candidate.sum())} candidates -> {out}')


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, help='prepared stage-2 input (lr.csv, entities.csv)')
    p.add_argument('--output', required=True)
    p.add_argument('--previous', help='earlier inventory directory; writes changes_vs_previous.csv')
    p.add_argument('--workers', type=int, default=2)
    return p


if __name__ == '__main__':
    run(parser().parse_args())
