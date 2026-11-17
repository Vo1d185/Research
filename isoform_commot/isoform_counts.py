"""UniProt isoform counts for the genes of the PLM LR benchmark (H-A: human, single-gene, all lengths).

Per gene: its reviewed UniProt entry (the benchmark mapping, dataset/sources/protein_map_human.csv) and the isoforms
listed in that entry's ALTERNATIVE PRODUCTS comment. A gene without that comment counts as 1 isoform.
  n_with_sequence  isoforms with a sequence in UniProt (status Displayed or Described): what PLM could score
  n_listed         every listed isoform (also 'Not described' and 'External')
  n_splicing       n_with_sequence if the comment's events include 'Alternative splicing', else 1
A pair has S' != 1 possible when ligand or receptor has >= 2 isoforms (the formula normalises over the isoform
combinations of one LR gene pair).

  fetch   download the comments (UniProt REST, generic client name) to data/lr_benchmark/isoforms/
  count   gene table + pair table + summary to reports/isoform_counts/
"""
import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import pandas as pd

from proteins import http_json
from research import check, launch_info

FIELDS = 'accession,id,gene_primary,length,cc_alternative_products'
WITH_SEQUENCE = {'Displayed', 'Described'}


def fetch(args):
    out = Path(args.output)
    check(not (out/'entries.json').exists(), f'{out}/entries.json exists')
    genes = pd.read_csv(args.protein_map, dtype=str)
    accessions = sorted(set(genes.accession))
    entries, releases = [], set()
    for i in range(0, len(accessions), 100):
        batch = accessions[i:i + 100]
        query = quote(' OR '.join(f'accession:{a}' for a in batch))
        url = f'https://rest.uniprot.org/uniprotkb/search?query={query}&fields={FIELDS}&format=json&size=500'
        data, headers = http_json(url)
        releases.add(headers.get('x-uniprot-release'))
        entries += data['results']
    found = {e['primaryAccession'] for e in entries}
    check(found == set(accessions), f'missing entries: {sorted(set(accessions) - found)[:10]}')
    out.mkdir(parents=True, exist_ok=True)
    (out/'entries.json').write_text(json.dumps(entries))
    (out/'fetch.json').write_text(json.dumps(dict(fetched_utc=datetime.now(timezone.utc).isoformat(),
                                                  uniprot_release=sorted(r for r in releases if r),
                                                  fields=FIELDS, accessions=len(accessions)), indent=2))


def isoform_row(entry):
    comments = [c for c in entry.get('comments', []) if c.get('commentType') == 'ALTERNATIVE PRODUCTS']
    check(len(comments) <= 1, f"{entry['primaryAccession']}: several ALTERNATIVE PRODUCTS comments")
    if not comments:
        return dict(n_listed=1, n_with_sequence=1, n_splicing=1, events='', isoform_ids='')
    c = comments[0]
    isoforms = c.get('isoforms', [])
    events = sorted(set(c.get('events', [])))
    seq = [i for i in isoforms if i.get('isoformSequenceStatus') in WITH_SEQUENCE]
    ids = [x for i in seq for x in i.get('isoformIds', [])]
    n_seq = max(len(seq), 1)
    return dict(n_listed=max(len(isoforms), 1), n_with_sequence=n_seq,
                n_splicing=n_seq if 'Alternative splicing' in events else 1,
                events=';'.join(events), isoform_ids=';'.join(ids))


def count(args):
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    entries = {e['primaryAccession']: e for e in json.loads((Path(args.isoforms)/'entries.json').read_text())}
    fetched = json.loads((Path(args.isoforms)/'fetch.json').read_text())
    genes = pd.read_csv(args.protein_map, dtype=str)
    table = pd.DataFrame([dict(gene=g, accession=a, **isoform_row(entries[a]))
                          for g, a in zip(genes.lr_symbol, genes.accession)])
    pairs = pd.read_csv(Path(args.dataset)/'positives.csv', dtype=str)
    check((pairs.single == 'True').all(), 'H-A is single-gene only')
    by_acc = table.drop_duplicates('accession').set_index('accession')
    rows = []
    for p in pairs.itertuples(index=False):
        lig, rec = by_acc.loc[p.ligand_accessions], by_acc.loc[p.receptor_accessions]
        row = dict(ligand=p.ligand, receptor=p.receptor, classes=p.classes, max_len=int(float(p.max_len)),
                   ligand_accession=p.ligand_accessions, receptor_accession=p.receptor_accessions)
        for k in ['n_with_sequence', 'n_listed', 'n_splicing']:
            row[f'ligand_{k}'], row[f'receptor_{k}'] = int(lig[k]), int(rec[k])
            row[f'combinations_{k}'] = int(lig[k]) * int(rec[k])
            row[f'multi_{k}'] = bool(int(lig[k]) >= 2 or int(rec[k]) >= 2)
        rows.append(row)
    pairs = pd.DataFrame(rows)
    summary = []
    for k in ['n_with_sequence', 'n_listed', 'n_splicing']:
        for subset, g in [('all', pairs), ('both <= 800 aa', pairs[pairs.max_len <= 800]),
                          ('any > 800 aa', pairs[pairs.max_len > 800])] + \
                         [(f'class {c}', pairs[pairs.classes == c]) for c in sorted(pairs.classes.unique())]:
            lig2, rec2 = g[f'ligand_{k}'] >= 2, g[f'receptor_{k}'] >= 2
            summary.append(dict(count=k, subset=subset, pairs=len(g), multi_either=int((lig2 | rec2).sum()),
                                fraction_multi_either=float((lig2 | rec2).mean()) if len(g) else float('nan'),
                                ligand_multi=int(lig2.sum()), receptor_multi=int(rec2.sum()), both_multi=int((lig2 & rec2).sum()),
                                median_combinations=float(g[f'combinations_{k}'].median()),
                                max_combinations=int(g[f'combinations_{k}'].max()),
                                total_combinations=int(g[f'combinations_{k}'].sum())))
    used = sorted(set(pairs.ligand_accession) | set(pairs.receptor_accession))
    gene_stats = table[table.accession.isin(used)].drop_duplicates('accession')
    out.mkdir(parents=True)
    table.to_csv(out/'genes.csv', index=False)
    pairs.to_csv(out/'pairs_H-A.csv', index=False)
    pd.DataFrame(summary).to_csv(out/'summary.csv', index=False)
    hist = gene_stats.n_with_sequence.value_counts().sort_index()
    (out/'manifest.json').write_text(json.dumps(dict(status='complete', launch=launch_info(), fetch=fetched,
        protein_map=args.protein_map, dataset=args.dataset, proteins_in_pairs=len(gene_stats),
        proteins_multi_with_sequence=int((gene_stats.n_with_sequence >= 2).sum()),
        isoforms_with_sequence_histogram={int(k): int(v) for k, v in hist.items()}, note=__doc__), indent=2, default=str))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fetch')
    f.add_argument('--protein-map', default='reports/plm_lr_benchmark/dataset/sources/protein_map_human.csv')
    f.add_argument('--output', default='data/lr_benchmark/isoforms')
    c = sub.add_parser('count')
    c.add_argument('--protein-map', default='reports/plm_lr_benchmark/dataset/sources/protein_map_human.csv')
    c.add_argument('--isoforms', default='data/lr_benchmark/isoforms')
    c.add_argument('--dataset', default='reports/plm_lr_benchmark/dataset/H-A')
    c.add_argument('--output', default='reports/isoform_counts')
    return p


if __name__ == '__main__':
    args = parser().parse_args()
    {'fetch': fetch, 'count': count}[args.command](args)
