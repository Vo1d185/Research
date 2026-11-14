"""Human-sequence PLM-interact diagnostics: pair list -> UniProt -> plm_input.csv (commot env).

Not part of the OT pipeline. Used to ask whether PLM-interact scores known
ligand-receptor pairs highly at all, independent of the mouse panel:

  * the mouse LR panel re-run with HUMAN orthologous sequences (species check);
  * classical human ligand-receptor positive controls, split by receptor length
    (<= 800 aa = inside the length range of the official test sets).

Each gene maps to exactly one reviewed human UniProt entry whose PRIMARY gene
name equals the symbol; anything else is refused unless an override with a
reason is given. Sequences are UniProt canonical, full length (precursors
included, e.g. preproinsulin, pro-TGF-beta1). Score the output with
`plm_score.py score --allow-long`, then run summarize.py.

Cross pairs (ligand x non-cognate receptor of the same set) are NOT negatives:
they are unannotated pairs used only to ask whether cognate pairs rank higher.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from urllib.parse import quote

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from proteins import (ESM_LETTERS, OFFICIAL_MAX_TOKENS, SEARCH, TRAINING_MAX_AA, UNIPROT,  # noqa: E402
                      digest, entry_fields, http_json, sha256_text, write_json)
from research import check  # noqa: E402

HUMAN_TAXON = 9606


def search_human(symbol, fetch=http_json):
    query = quote(f'(gene_exact:{symbol}) AND (organism_id:{HUMAN_TAXON}) AND (reviewed:true)')
    found, headers = fetch(SEARCH.format(query=query))
    hits = [(e, f) for e in found.get('results', []) for f in [entry_fields(e)]
            if f['taxon_id'] == HUMAN_TAXON and f['reviewed']
            and str(f['gene_primary']).upper() == symbol.upper()]
    return hits, headers


def resolve_human(symbol, override=None, fetch=http_json):
    if override:
        entry, headers = fetch(f'{UNIPROT}/{quote(override)}.json')
        f = entry_fields(entry)
        check(f['taxon_id'] == HUMAN_TAXON and f['reviewed'], f'{override} is not a reviewed human entry')
        check(str(f['gene_primary']).upper() == symbol.upper(),
              f'Override {override} has primary gene {f["gene_primary"]}, not {symbol}')
        return f, headers, 'manual_override'
    hits, headers = search_human(symbol, fetch)
    check(len(hits) == 1, f'{symbol}: {len(hits)} reviewed human entries with this primary gene name '
          f'{[h[1]["accession"] for h in hits]}; add a row to the overrides file (gene,accession,reason)')
    return hits[0][1], headers, 'unique_reviewed_primary_gene'


def build(pairs, overrides=None, known=None, cross=True, fetch=http_json):
    """Return (protein_map, jobs) DataFrames; jobs include both input orders."""
    pairs = pairs[['ligand', 'receptor']].astype(str)
    check(len(pairs) and not pairs.duplicated().any(), 'Empty or duplicate pair list')
    overrides = overrides or {}
    known = known or {}
    genes = sorted(set(pairs.ligand) | set(pairs.receptor))
    unknown = sorted(set(overrides) - set(genes))
    check(not unknown, f'Overrides for genes not in the pair list: {unknown}')
    rows, seqs, errors = [], {}, []
    for gene in genes:
        try:
            f, headers, how = resolve_human(gene, overrides.get(gene), fetch)
        except ValueError as exc:  # report every ambiguous gene at once
            errors.append(str(exc))
            continue
        bad = sorted(set(f['sequence']) - ESM_LETTERS)
        check(f['sequence'] and not bad, f'{gene}: empty sequence or letters outside ESM vocabulary {bad}')
        seqs[gene] = f['sequence']
        rows.append(dict(gene=gene, accession=f['accession'], uniprot_id=f['uniprot_id'],
            gene_primary=f['gene_primary'], protein_name=f['protein_name'], taxon_id=f['taxon_id'],
            isoform='canonical', entry_version=f['entry_version'], sequence_version=f['sequence_version'],
            length=len(f['sequence']), sequence_sha256=sha256_text(f['sequence']),
            uniprot_release=headers.get('x-uniprot-release'), resolution=how,
            retrieved_utc=datetime.now(timezone.utc).isoformat(timespec='seconds')))
    check(not errors, ' | '.join(errors))
    pmap = pd.DataFrame(rows)
    shared = pmap[pmap.accession.duplicated(keep=False)]
    check(shared.empty, f'Different genes mapped to the same accession: {shared.gene.tolist()}')
    acc = dict(zip(pmap.gene, pmap.accession))
    sets = [('cognate', r.ligand, r.receptor) for r in pairs.itertuples(index=False)]
    if cross:
        cognate = set(zip(pairs.ligand, pairs.receptor))
        sets += [('cross', l, r) for l in sorted(set(pairs.ligand)) for r in sorted(set(pairs.receptor))
                 if (l, r) not in cognate and l != r]
    jobs = []
    for kind, ligand, receptor in sets:
        group = 'receptor_le_800' if len(seqs[receptor]) <= TRAINING_MAX_AA else 'receptor_gt_800'
        for order, (q, t) in [('LR', (ligand, receptor)), ('RL', (receptor, ligand))]:
            jobs.append(dict(pair_id=f'{ligand}__{receptor}__{order}', set=kind, ligand=ligand,
                receptor=receptor, order=order, receptor_group=group,
                ligand_len=len(seqs[ligand]), receptor_len=len(seqs[receptor]),
                known_cross_reactivity=known.get((ligand, receptor)) if kind == 'cross' else None,
                query_symbol=q, text_symbol=t, query_accession=acc[q], text_accession=acc[t],
                n_tokens=len(seqs[q]) + len(seqs[t]) + 3,
                beyond_training_length=max(len(seqs[q]), len(seqs[t])) > TRAINING_MAX_AA,
                query=seqs[q], text=seqs[t]))
    jobs = pd.DataFrame(jobs)
    check(jobs.pair_id.is_unique, 'Duplicate job IDs')
    return pmap, jobs, seqs


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--pairs', required=True, help='CSV ligand,receptor')
    p.add_argument('--overrides', help='CSV gene,accession,reason')
    p.add_argument('--known-cross', help='CSV ligand,receptor,note: documented cross-reactive pairs')
    p.add_argument('--no-cross', action='store_true')
    p.add_argument('--output', required=True)
    args = p.parse_args()
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    overrides = {}
    if args.overrides:
        o = pd.read_csv(args.overrides, dtype=str)
        check({'gene', 'accession', 'reason'} <= set(o) and o.reason.notna().all() and o.gene.is_unique,
              'Overrides need gene,accession,reason (one row per gene, reason required)')
        overrides = dict(zip(o.gene, o.accession))
    known = {}
    if args.known_cross:
        k = pd.read_csv(args.known_cross, dtype=str)
        known = {(r.ligand, r.receptor): r.note for r in k.itertuples(index=False)}
    pairs = pd.read_csv(args.pairs, dtype=str)
    pmap, jobs, seqs = build(pairs, overrides, known, not args.no_cross)
    out.mkdir(parents=True)
    pmap.to_csv(out/'protein_map.csv', index=False)
    with (out/'proteins.fasta').open('w') as fasta:
        for row in pmap.itertuples(index=False):
            fasta.write(f'>{row.accession}|{row.gene}|human\n{seqs[row.gene]}\n')
    jobs.drop(columns=['query', 'text']).to_csv(out/'jobs.csv', index=False)
    jobs[['pair_id', 'query', 'text']].to_csv(out/'plm_input.csv', index=False)
    cognate = jobs[(jobs.set == 'cognate') & (jobs.order == 'LR')]
    write_json(out/'jobs.json', dict(status='complete', species='Homo sapiens (taxon 9606)',
        sequence='UniProt canonical, full length (precursor incl. signal peptide/propeptides)',
        pairs_file=str(Path(args.pairs)), pairs_sha256=digest(args.pairs),
        overrides_sha256=digest(args.overrides) if args.overrides else None,
        known_cross_sha256=digest(args.known_cross) if args.known_cross else None,
        uniprot_releases=sorted(set(pmap.uniprot_release.dropna())),
        n_jobs=len(jobs), n_cognate_pairs=len(cognate), n_cross_pairs=int((jobs.set == 'cross').sum() // 2),
        cognate_by_receptor_group=cognate.receptor_group.value_counts().to_dict(),
        grouping=f'receptor UniProt canonical length <= {TRAINING_MAX_AA} aa vs > {TRAINING_MAX_AA} aa',
        max_tokens=int(jobs.n_tokens.max()),
        jobs_over_official_1603_tokens=jobs.loc[jobs.n_tokens > OFFICIAL_MAX_TOKENS, 'pair_id'].tolist(),
        plm_input_sha256=digest(out/'plm_input.csv'), fasta_sha256=digest(out/'proteins.fasta'),
        note='cross pairs are unannotated, not proven non-interacting'))
    print(f'{len(pmap)} proteins, {len(jobs)} jobs -> {out}/plm_input.csv; max tokens {jobs.n_tokens.max()}')
    print(cognate[['ligand', 'receptor', 'ligand_len', 'receptor_len', 'receptor_group']].to_string(index=False))


if __name__ == '__main__':
    main()
