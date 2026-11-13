"""Gene -> mouse protein -> PLM-interact jobs -> ppi.csv (runs in the commot env).

  map     LR symbols of a prepared input -> original data symbols -> exactly one
          reviewed (Swiss-Prot) MOUSE UniProt entry per gene; canonical sequence.
  pairs   PLM-interact jobs: every cognate LR pair in both orders, plus panel
          cross pairs (ligand x receptor not in lr.csv) used ONLY to check the scores.
  import  PLM scores -> ppi.csv for the cognate pairs (mean of both orders) and a
          score report. Never imputes a missing score.

PLM-interact scores are model predictions from sequence, not measured affinities.
The released checkpoints were trained on human PPIs; mouse use is cross-species.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

from research import check

MOUSE_TAXON = 10090
ESM_LETTERS = set('ACDEFGHIKLMNPQRSTVWYXBUZO')
TRAINING_MAX_AA = 800          # official test sets contain 50-800 aa proteins
OFFICIAL_MAX_TOKENS = 1603     # 800 + 800 + 3 special tokens
UNIPROT = 'https://rest.uniprot.org/uniprotkb'
SEARCH = UNIPROT + '/search?query={query}&format=json&size=25'


def sha256_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False))


def http_json(url):
    """GET JSON from UniProt. No personal data is sent; generic client name only."""
    request = Request(url, headers={'Accept': 'application/json',
                                    'User-Agent': 'isoform-commot-research/1.0'})
    with urlopen(request, timeout=60) as response:
        headers = {k.lower(): v for k, v in response.headers.items()}
        return json.loads(response.read().decode()), headers


# ---------------------------------------------------------------- map
def data_symbols(input_dir, h5ad=None, gene_map=None, symbol_table=None):
    """LR symbol (as in lr.csv) -> symbol in the expression data (e.g. APOE -> Apoe)."""
    lr = pd.read_csv(Path(input_dir)/'lr.csv', dtype=str)
    wanted = sorted(set(lr.ligand) | set(lr.receptor))
    if symbol_table:
        table = pd.read_csv(symbol_table, dtype=str)
        check({'lr_symbol', 'data_symbol'} <= set(table), 'symbol table needs lr_symbol,data_symbol')
        check(not table.lr_symbol.duplicated().any(), 'Duplicate lr_symbol in symbol table')
        mapping = dict(zip(table.lr_symbol, table.data_symbol))
        source = f'symbol_table sha256={digest(symbol_table)}'
    else:
        import stage1
        prep = json.loads((Path(input_dir)/'preparation.json').read_text())
        check(prep.get('status') == 'complete', 'Need a complete prepared input')
        h5ad = h5ad or prep['h5ad']
        if prep.get('gene_map_sha256'):
            check(gene_map is not None and digest(gene_map) == prep['gene_map_sha256'],
                  'This input used --gene-map; pass the same file with --gene-map')
        data = stage1.open_data(h5ad)
        try:
            _, genes = stage1.matrix_and_genes(data, prep['matrix'], prep.get('gene_key'))
        finally:
            data.file.close()
        mapped = stage1.map_genes(genes, prep['symbol_policy'], gene_map)
        mapping = {m: g for g, m in zip(genes, mapped) if m is not None}
        source = f'h5ad var ({prep["matrix"]}, gene_key={prep.get("gene_key")}, policy={prep["symbol_policy"]})'
    missing = [g for g in wanted if g not in mapping]
    check(not missing, f'No data symbol for {missing}')
    return {g: mapping[g] for g in wanted}, source


def entry_fields(entry):
    genes = entry.get('genes') or []
    primary = genes[0].get('geneName', {}).get('value') if genes else None
    audit = entry.get('entryAudit', {})
    sequence = entry.get('sequence', {}).get('value', '')
    return dict(accession=entry['primaryAccession'], uniprot_id=entry.get('uniProtkbId'),
        gene_primary=primary, taxon_id=entry.get('organism', {}).get('taxonId'),
        reviewed='Swiss-Prot' in str(entry.get('entryType', '')),
        protein_name=entry.get('proteinDescription', {}).get('recommendedName', {})
                          .get('fullName', {}).get('value'),
        entry_version=audit.get('entryVersion'), sequence_version=audit.get('sequenceVersion'),
        sequence=sequence)


def resolve(symbol, override=None, fetch=http_json):
    """Exactly one reviewed mouse entry whose PRIMARY gene name is the symbol."""
    if override:
        entry, headers = fetch(f'{UNIPROT}/{quote(override)}.json')
        fields = entry_fields(entry)
        check(fields['taxon_id'] == MOUSE_TAXON, f'Override {override} is not a mouse entry')
        return entry, fields, headers, 'manual_override'
    query = quote(f'(gene_exact:{symbol}) AND (organism_id:{MOUSE_TAXON}) AND (reviewed:true)')
    found, headers = fetch(SEARCH.format(query=query))
    hits = []
    for entry in found.get('results', []):
        fields = entry_fields(entry)
        if (fields['taxon_id'] == MOUSE_TAXON and fields['reviewed'] and
                str(fields['gene_primary']).lower() == symbol.lower()):
            hits.append((entry, fields))
    check(len(hits) == 1, f'{symbol}: {len(hits)} reviewed mouse entries with this primary gene name '
          f'{[h[1]["accession"] for h in hits]}; resolve manually with --override (lr_symbol,accession,reason)')
    return hits[0][0], hits[0][1], headers, 'unique_reviewed_primary_gene'


def map_proteins(args, fetch=http_json):
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    symbols, symbol_source = data_symbols(args.input, args.h5ad, args.gene_map, args.symbol_table)
    overrides = {}
    if args.override:
        table = pd.read_csv(args.override, dtype=str)
        check({'lr_symbol', 'accession', 'reason'} <= set(table) and table.reason.notna().all(),
              'Override needs lr_symbol,accession,reason (reason is required)')
        overrides = dict(zip(table.lr_symbol, table.accession))
    rows, raw, releases = [], {}, set()
    for lr_symbol, data_symbol in symbols.items():
        entry, f, headers, how = resolve(data_symbol, overrides.get(lr_symbol), fetch)
        bad = sorted(set(f['sequence']) - ESM_LETTERS)
        check(f['sequence'] and not bad, f'{lr_symbol}: empty sequence or letters outside ESM vocabulary {bad}')
        release = headers.get('x-uniprot-release')
        releases.add(release)
        raw[f['accession']] = entry
        rows.append(dict(lr_symbol=lr_symbol, data_symbol=data_symbol, accession=f['accession'],
            uniprot_id=f['uniprot_id'], gene_primary=f['gene_primary'], protein_name=f['protein_name'],
            taxon_id=f['taxon_id'], isoform='canonical', entry_version=f['entry_version'],
            sequence_version=f['sequence_version'], length=len(f['sequence']),
            sequence_sha256=sha256_text(f['sequence']), uniprot_release=release,
            uniprot_release_date=headers.get('x-uniprot-release-date'),
            retrieved_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'),
            resolution=how, longer_than_training_max=len(f['sequence']) > TRAINING_MAX_AA,
            sequence=f['sequence']))
    table = pd.DataFrame(rows)
    shared = table[table.accession.duplicated(keep=False)]
    check(shared.empty, f'Different genes mapped to the same accession: {shared.lr_symbol.tolist()}')
    out.mkdir(parents=True)
    (out/'uniprot').mkdir()
    for accession, entry in raw.items():
        write_json(out/'uniprot'/f'{accession}.json', entry)
    table.drop(columns='sequence').to_csv(out/'protein_map.csv', index=False)
    with (out/'proteins.fasta').open('w') as fasta:
        for row in table.itertuples(index=False):
            fasta.write(f'>{row.accession}|{row.lr_symbol}|{row.data_symbol}\n{row.sequence}\n')
    write_json(out/'proteins.json', dict(status='complete', input=str(Path(args.input).resolve()),
        lr_sha256=digest(Path(args.input)/'lr.csv'), symbol_source=symbol_source,
        species='Mus musculus (taxon 10090)', sequence='UniProt canonical, full length incl. signal peptide',
        uniprot_releases=sorted(r for r in releases if r), n_proteins=len(table),
        longer_than_training_max=table.loc[table.longer_than_training_max, 'lr_symbol'].tolist(),
        rule='exactly one reviewed mouse entry whose primary gene name equals the data symbol; '
             'otherwise refused unless a documented override is given',
        fasta_sha256=digest(out/'proteins.fasta')))
    print(f'Mapped {len(table)} genes; review {out}/protein_map.csv before scoring')


def read_fasta(path):
    seqs, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split('|')[1]
            seqs[name] = ''
        elif name:
            seqs[name] += line.strip()
    return seqs


# ---------------------------------------------------------------- pairs
def make_pairs(args):
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    proteins = Path(args.proteins)
    meta = json.loads((proteins/'proteins.json').read_text())
    check(meta['status'] == 'complete' and digest(proteins/'proteins.fasta') == meta['fasta_sha256'],
          'Protein map incomplete or changed')
    check(digest(Path(args.input)/'lr.csv') == meta['lr_sha256'], 'lr.csv differs from the mapped panel')
    seqs = read_fasta(proteins/'proteins.fasta')
    pmap = pd.read_csv(proteins/'protein_map.csv', dtype=str).set_index('lr_symbol')
    lr = pd.read_csv(Path(args.input)/'lr.csv', dtype=str)
    database = set()
    if args.lr_dir:
        import stage1
        full = stage1.load_lr(args.lr_dir)
        database = set(zip(full.ligand, full.receptor))
    sets = [('cognate', row.ligand, row.receptor, row.type) for row in lr.itertuples(index=False)]
    if not args.no_cross:
        cognate = set(zip(lr.ligand, lr.receptor))
        for ligand in sorted(set(lr.ligand)):
            for receptor in sorted(set(lr.receptor)):
                if (ligand, receptor) not in cognate and ligand != receptor:
                    sets.append(('panel_cross', ligand, receptor, None))
    jobs = []
    for kind, ligand, receptor, lr_type in sets:
        for order, (q, t) in [('LR', (ligand, receptor)), ('RL', (receptor, ligand))]:
            jobs.append(dict(pair_id=f'{ligand}__{receptor}__{order}', set=kind, ligand=ligand,
                receptor=receptor, type=lr_type, order=order,
                in_lr_database=(ligand, receptor) in database if args.lr_dir else None,
                query_symbol=q, text_symbol=t, query_accession=pmap.loc[q, 'accession'],
                text_accession=pmap.loc[t, 'accession'], query_len=len(seqs[q]), text_len=len(seqs[t]),
                n_tokens=len(seqs[q]) + len(seqs[t]) + 3,
                beyond_training_length=max(len(seqs[q]), len(seqs[t])) > TRAINING_MAX_AA,
                query=seqs[q], text=seqs[t]))
    jobs = pd.DataFrame(jobs)
    check(jobs.pair_id.is_unique, 'Duplicate job IDs')
    out.mkdir(parents=True)
    jobs.to_csv(out/'jobs.csv', index=False)
    jobs[['pair_id', 'query', 'text']].to_csv(out/'plm_input.csv', index=False)
    write_json(out/'jobs.json', dict(status='complete', proteins_fasta_sha256=meta['fasta_sha256'],
        lr_sha256=meta['lr_sha256'], n_jobs=len(jobs),
        n_cognate_pairs=int((jobs.set == 'cognate').sum() // 2),
        n_cross_pairs=int((jobs.set == 'panel_cross').sum() // 2),
        cross_pairs_in_lr_database=int(jobs.in_lr_database.fillna(False).astype(bool).sum() // 2),
        max_tokens=int(jobs.n_tokens.max()),
        jobs_over_official_1603_tokens=jobs.loc[jobs.n_tokens > OFFICIAL_MAX_TOKENS, 'pair_id'].tolist(),
        jobs_beyond_training_length=int(jobs.beyond_training_length.sum()),
        plm_input_sha256=digest(out/'plm_input.csv'),
        note='panel_cross pairs are NOT negatives: unannotated, not proven non-interacting; '
             'they only check whether cognate pairs score higher.'))
    print(f'{len(jobs)} jobs -> {out}/plm_input.csv; max tokens {jobs.n_tokens.max()}')


# ---------------------------------------------------------------- import
def auroc(positive, negative):
    """Mann-Whitney AUROC with ties counted as 1/2; None when a group is empty."""
    positive, negative = np.asarray(positive, float), np.asarray(negative, float)
    if not len(positive) or not len(negative):
        return None
    greater = (positive[:, None] > negative[None, :]).sum()
    ties = (positive[:, None] == negative[None, :]).sum()
    return float((greater + 0.5*ties) / (len(positive)*len(negative)))


def import_scores(args):
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    jobs_dir = Path(args.jobs)
    jobs_meta = json.loads((jobs_dir/'jobs.json').read_text())
    jobs = pd.read_csv(jobs_dir/'jobs.csv', dtype={'ligand': str, 'receptor': str})
    check(digest(jobs_dir/'plm_input.csv') == jobs_meta['plm_input_sha256'], 'Jobs changed after export')
    scores = pd.read_csv(args.scores, dtype={'pair_id': str})
    meta = json.loads(Path(args.scores_meta or str(args.scores) + '.meta.json').read_text())
    check(meta.get('status') == 'complete' and meta.get('input_sha256') == jobs_meta['plm_input_sha256'],
          'Scores were not produced from these jobs (input hash mismatch)')
    check(scores.pair_id.is_unique and set(scores.pair_id) == set(jobs.pair_id),
          'Scores must cover every job exactly once; no imputation')
    s = scores.score.to_numpy(float)
    check(np.isfinite(s).all() and ((s >= 0) & (s <= 1)).all(), 'Scores must be in [0,1]')
    merged = jobs.merge(scores[['pair_id', 'score', 'n_tokens']], on='pair_id', suffixes=('', '_model'))
    check((merged.n_tokens == merged.n_tokens_model).all(), 'Token counts differ: sequences were truncated/altered')
    wide = merged.pivot_table(index=['set', 'ligand', 'receptor'], columns='order', values='score').reset_index()
    wide['score'] = (wide.LR + wide.RL) / 2      # pre-specified: mean of both orders
    wide['abs_order_diff'] = (wide.LR - wide.RL).abs()
    db = merged.groupby(['ligand', 'receptor']).in_lr_database.first()
    wide['in_lr_database'] = [db.get((l, r)) for l, r in zip(wide.ligand, wide.receptor)]
    lr = pd.read_csv(Path(args.input)/'lr.csv', dtype=str)
    check(digest(Path(args.input)/'lr.csv') == jobs_meta['lr_sha256'], 'lr.csv differs from the scored panel')
    cognate = lr.merge(wide[wide.set == 'cognate'], on=['ligand', 'receptor'], how='left', validate='one_to_one')
    check(cognate.score.notna().all(), 'Missing cognate score')
    cross = wide[wide.set == 'panel_cross']
    unannotated = cross[cross.in_lr_database != True]  # noqa: E712 (None means database not checked)
    out.mkdir(parents=True)
    cognate.rename(columns={'LR': 'score_LR', 'RL': 'score_RL'})[
        ['ligand', 'receptor', 'type', 'score', 'score_LR', 'score_RL', 'abs_order_diff']].to_csv(out/'ppi.csv', index=False)
    wide.rename(columns={'LR': 'score_LR', 'RL': 'score_RL'}).to_csv(out/'scores_by_pair.csv', index=False)
    provenance = (f'PLM-interact {meta.get("checkpoint_name")}@{meta.get("checkpoint_revision")} '
                  f'(pytorch_model.bin sha256 {str(meta.get("checkpoint_sha256"))[:16]}); '
                  f'mouse UniProt canonical sequences; mean of both input orders; no truncation; '
                  f'model prediction, not measured affinity')
    report = dict(status='complete', score_provenance=provenance, model=meta,
        jobs_sha256=jobs_meta['plm_input_sha256'], ppi_sha256=digest(out/'ppi.csv'),
        cognate=dict(n=int(len(cognate)), min=float(cognate.score.min()), max=float(cognate.score.max()),
                     mean=float(cognate.score.mean()),
                     by_type=cognate.groupby('type').score.agg(['min', 'mean', 'max']).to_dict('index')),
        order_asymmetry=dict(max=float(wide.abs_order_diff.max()), mean=float(wide.abs_order_diff.mean())),
        cognate_vs_unannotated_cross_auroc=auroc(cognate.score, unannotated.score),
        n_unannotated_cross=int(len(unannotated)),
        jobs_beyond_training_length=jobs_meta['jobs_beyond_training_length'],
        interpretation=['Scores near 1 for all cognate pairs mean the PPI cost term is ~0 (saturated prior).',
                        'AUROC uses unannotated panel pairs, which are not proven negatives.',
                        'Checkpoints are trained on human PPIs; mouse use is cross-species.'])
    write_json(out/'score_report.json', report)
    print(json.dumps({k: report[k] for k in ['cognate', 'order_asymmetry',
                                            'cognate_vs_unannotated_cross_auroc']}, indent=2))
    print(f'ppi.csv -> {out}; attach with: python stage1.py attach-ppi --ppi-dir {out} ...')


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    m = sub.add_parser('map')
    m.add_argument('--input', required=True, help='prepared input directory (stage1_input_01)')
    m.add_argument('--output', required=True)
    m.add_argument('--h5ad', help='defaults to the path in preparation.json')
    m.add_argument('--gene-map', help='same --gene-map file used by prepare, if any')
    m.add_argument('--symbol-table', help='CSV lr_symbol,data_symbol instead of reading the h5ad')
    m.add_argument('--override', help='CSV lr_symbol,accession,reason for manually resolved genes')
    q = sub.add_parser('pairs')
    q.add_argument('--input', required=True); q.add_argument('--proteins', required=True)
    q.add_argument('--output', required=True)
    q.add_argument('--lr-dir', help='pinned CCCvelo LR tables; marks cross pairs that are annotated')
    q.add_argument('--no-cross', action='store_true', help='score cognate pairs only')
    i = sub.add_parser('import')
    i.add_argument('--input', required=True); i.add_argument('--jobs', required=True)
    i.add_argument('--scores', required=True); i.add_argument('--scores-meta')
    i.add_argument('--output', required=True)
    return p


def main():
    args = parser().parse_args()
    {'map': map_proteins, 'pairs': make_pairs, 'import': import_scores}[args.command](args)


if __name__ == '__main__':
    main()
