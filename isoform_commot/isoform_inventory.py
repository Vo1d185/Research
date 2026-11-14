"""Isoform inventory of the LR panel genes (UniProt, mouse). Counts and flags only; no analysis.

For every gene of a prepared stage-2 input (ligands and receptors; complexes split into subunits):
reviewed mouse UniProt entries whose primary gene name is the gene (proteins.search), their isoforms
(ALTERNATIVE PRODUCTS: canonical + alternative products) and 'Alternative sequence' features
(variable regions, canonical coordinates). Each variable region is flagged when it overlaps

  * the extracellular region: 'Topological domain' features described as Extracellular; for an entry
    without topological domains whose subcellular location includes Secreted, the chain after the
    signal peptide (the whole chain if no signal peptide is annotated);
  * an annotated domain ('Domain' features);
  * an annotated binding region ('Region' features whose description contains 'bind' or 'interact')
    or a 'Binding site' feature.

A gene is a candidate when one of its entries has >= 2 isoforms and a variable region with any flag.
Genes with several reviewed entries (e.g. Nrxn1 alpha/beta, separate entries) are listed as such;
Nrxn1 alpha/beta is the first row of candidates.csv as the known example.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time
from urllib.error import HTTPError

import pandas as pd

import proteins
from research import check, launch_info

KNOWN = dict(gene='Nrxn1', note='known example: Nrxn1 alpha (Q9CS84) and beta (P0DI97) are separate reviewed '
             'entries (alternative promoters); see HANDOFF.md (beta used in stage 1)')


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


def describe(gene, entry, fields):
    length = len(fields['sequence'])
    ec, basis = extracellular(entry, length)
    domains = labelled(entry, 'Domain')
    binding = labelled(entry, 'Region', lambda d: 'bind' in d.lower() or 'interact' in d.lower()) + \
        labelled(entry, 'Binding site')
    variants = []
    for start, end, description in labelled(entry, 'Alternative sequence'):
        region = (start, end)
        hit_domains = overlaps(region, [(s, e) for s, e, _ in domains])
        hit_binding = overlaps(region, [(s, e) for s, e, _ in binding])
        variants.append(dict(gene=gene, accession=fields['accession'], start=start, end=end,
            description=description, overlaps_extracellular=bool(overlaps(region, ec)),
            overlaps_domain=bool(hit_domains), overlaps_binding=bool(hit_binding),
            domains=';'.join(sorted({d for s, e, d in domains if overlaps(region, [(s, e)])})),
            binding=';'.join(sorted({d for s, e, d in binding if overlaps(region, [(s, e)])}))))
    iso = isoforms(entry)
    return dict(gene=gene, accession=fields['accession'], uniprot_id=fields['uniprot_id'],
                protein_name=fields['protein_name'], length=length, isoforms=max(len(iso), 1),
                isoform_names=';'.join(str(i['name']) for i in iso), variable_regions=len(variants),
                extracellular_basis=basis, extracellular_ranges=';'.join(f'{s}-{e}' for s, e in ec),
                entry_version=fields['entry_version'], sequence_sha256=proteins.sha256_text(fields['sequence'])), variants


def run(args, fetch=None):
    launch = launch_info()
    fetch = fetch or polite(proteins.http_json)
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    genes, lr = panel_genes(args.input)
    names = sorted(genes)

    def one(gene):
        hits, headers = proteins.search(gene, fetch)
        return gene, hits, headers
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(one, names))
    entries, variants, raw, releases = [], [], {}, set()
    for gene, hits, headers in results:
        releases.add(headers.get('x-uniprot-release'))
        if not hits:
            entries.append(dict(gene=gene, accession=None, protein_name='no reviewed mouse entry with this primary gene name'))
        for entry, fields in hits:
            row, var = describe(gene, entry, fields)
            entries.append(row)
            variants.extend(var)
            raw[fields['accession']] = proteins.sha256_text(json.dumps(entry, sort_keys=True))
    entries, variants = pd.DataFrame(entries), pd.DataFrame(variants)
    if variants.empty:
        variants = pd.DataFrame(columns=['gene', 'accession', 'overlaps_extracellular', 'overlaps_domain',
                                         'overlaps_binding'])
    flags = variants.groupby(['gene', 'accession'])[['overlaps_extracellular', 'overlaps_domain',
                                                      'overlaps_binding']].any().reset_index()
    entries = entries.merge(flags, on=['gene', 'accession'], how='left')
    for c in ['overlaps_extracellular', 'overlaps_domain', 'overlaps_binding']:
        entries[c] = entries[c].astype('boolean').fillna(False).astype(bool)
    per_gene = entries.groupby('gene').agg(
        reviewed_entries=('accession', 'count'), accessions=('accession', lambda s: ';'.join(s.dropna())),
        max_isoforms=('isoforms', 'max'), variable_regions=('variable_regions', 'sum'),
        any_extracellular=('overlaps_extracellular', 'any'), any_domain=('overlaps_domain', 'any'),
        any_binding=('overlaps_binding', 'any')).reset_index()
    per_gene['candidate'] = (per_gene.max_isoforms >= 2) & (per_gene.any_extracellular | per_gene.any_domain |
                                                            per_gene.any_binding)
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
    out.mkdir(parents=True)
    per_gene.to_csv(out/'genes.csv', index=False)
    entries.to_csv(out/'entries.csv', index=False)
    variants.to_csv(out/'variable_regions.csv', index=False)
    candidates.to_csv(out/'candidates.csv', index=False)
    prep = json.loads((Path(args.input)/'preparation.json').read_text())
    (out/'inventory.json').write_text(json.dumps(dict(status='complete', launch=launch,
        created=datetime.now(timezone.utc).isoformat(timespec='seconds'), input=str(args.input),
        sample=prep.get('sample'), input_sha256={f: hashlib.sha256((Path(args.input)/f).read_bytes()).hexdigest()
                                                 for f in ['lr.csv', 'entities.csv']},
        uniprot_releases=sorted(r for r in releases if r), query='proteins.search: (gene_exact:G) AND '
        '(organism_id:10090) AND (reviewed:true), primary gene name must equal G',
        genes=len(per_gene), genes_with_entry=int((per_gene.reviewed_entries > 0).sum()),
        genes_with_several_entries=sorted(per_gene[per_gene.reviewed_entries > 1].gene),
        candidates=int(per_gene.candidate.sum()), entry_json_sha256=raw), indent=2, default=str))
    print(f'{len(per_gene)} genes, {int(per_gene.candidate.sum())} candidates -> {out}')


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True, help='prepared stage-2 input (lr.csv, entities.csv)')
    p.add_argument('--output', required=True)
    p.add_argument('--workers', type=int, default=2)
    return p


if __name__ == '__main__':
    run(parser().parse_args())
