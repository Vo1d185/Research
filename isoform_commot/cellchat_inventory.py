"""CellChatDB (mouse) inventory on the whole tissue — a REPORT, not a panel selection.

  fetch   download the pinned CellChatDB.mouse.rda (CellChat v2, GPL-3, not
          redistributed), verify SHA256, convert interaction/complex tables to CSV
  report  for Secreted Signaling and Cell-Cell Contact separately, count what
          remains after each step, and list CellChat vs CCCvelo type conflicts:
            0  all CellChat pairs
            a  ligand and receptor are single genes (complexes counted separately)
            b  both genes detected (count > 0) in >= 5% of ALL spots of the h5ad
            c  receptor has a scSeqComm receptor->TF link AND, with the top-10 TF rule
               of evaluate_targets.py, >= 5 target genes detected in >= 5% of spots

Target exclusion in (c): every single-gene ligand/receptor of the pairs that pass (b)
in EITHER class is removed from the target lists. Any panel drawn from that pool
excludes a subset of these genes, so a receptor that passes here stays evaluable
for every such panel (conservative). Excluding only the pair's own ligand and
receptor is reported as a sensitivity count.

CellChat vs CCCvelo: CCCvelo tables use upper-case symbols; CellChat mouse symbols
are upper-cased for matching. This is symbol matching, not ortholog mapping.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from research import check

CELLCHAT_COMMIT = '75253cd0c9e68410e6e721a6d3a0419a1d7e358f'
CELLCHAT_URL = (f'https://raw.githubusercontent.com/jinworks/CellChat/{CELLCHAT_COMMIT}'
                '/data/CellChatDB.mouse.rda')
CELLCHAT_SHA256 = 'dd59b19646e546b3ed41f394ad2cf64cf44f3d3fcfbb58ac176f0fe1322f9413'
CELLCHAT_FILE_COMMIT = '592ba49cf4'  # last commit touching the file (2025-04-16), via GitHub API
CLASSES = ['Secreted Signaling', 'Cell-Cell Contact']
CCC_TYPE = {'Secreted Signaling': 'diffusion', 'Cell-Cell Contact': 'contact'}
MIN_DETECTED, TOP_TF, MIN_TARGETS = 0.05, 10, 5


def fetch(output):
    import rdata
    from stage1 import digest, download, write_json
    out = Path(output)
    check(not (out/'db.json').exists(), f'{out} already has a converted database')
    raw = out/'raw'/'CellChatDB.mouse.rda'
    download(CELLCHAT_URL, raw, expected=CELLCHAT_SHA256)
    obj = next(iter(rdata.read_rda(raw).values()))
    interaction = obj['interaction'].copy()
    interaction.columns = [str(c) for c in interaction.columns]
    complex_ = obj['complex'].copy()
    complex_.index = complex_.index.astype(str)
    complex_.columns = [str(c) for c in complex_.columns]
    genes = obj['geneInfo'][['Symbol']].astype(str)
    interaction.to_csv(out/'interaction.csv', index=False)
    complex_.to_csv(out/'complex.csv', index_label='complex')
    genes.to_csv(out/'gene_symbols.csv', index=False)
    write_json(out/'db.json', dict(status='complete', source='CellChat v2 (jinworks/CellChat)',
        commit=CELLCHAT_COMMIT, file_last_changed_commit=CELLCHAT_FILE_COMMIT, url=CELLCHAT_URL,
        raw_sha256=CELLCHAT_SHA256, sha256_source='calculated locally at first download, not author-published',
        license='GPL-3 (downloaded, not redistributed)',
        retrieved_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'),
        n_interactions=len(interaction), annotation_counts=interaction.annotation.value_counts().to_dict(),
        version_counts=interaction['version'].value_counts().to_dict(),
        converted_sha256={f: digest(out/f) for f in ['interaction.csv', 'complex.csv', 'gene_symbols.csv']}))
    print(f'CellChatDB mouse: {len(interaction)} interactions -> {out}')


def load_cellchat(folder):
    from stage1 import digest
    folder = Path(folder)
    meta = json.loads((folder/'db.json').read_text())
    check(meta.get('status') == 'complete' and meta['raw_sha256'] == CELLCHAT_SHA256, 'Run fetch first')
    for name, sha in meta['converted_sha256'].items():
        check(digest(folder/name) == sha, f'{name} changed after conversion')
    it = pd.read_csv(folder/'interaction.csv', dtype=str)
    complexes = set(pd.read_csv(folder/'complex.csv', dtype=str).complex)
    symbols = set(pd.read_csv(folder/'gene_symbols.csv', dtype=str).Symbol)
    for side in ['ligand', 'receptor']:
        unknown = ~it[side].isin(complexes) & ~it[side].isin(symbols)
        check(not unknown.any(), f'{side} entities neither complex nor gene: {it.loc[unknown, side].head().tolist()}')
    check(not it.duplicated(['ligand', 'receptor']).any(), 'Duplicate CellChat ligand-receptor pairs')
    it['ligand_complex'] = it.ligand.isin(complexes)
    it['receptor_complex'] = it.receptor.isin(complexes)
    return it, meta


def detected_fractions(h5ad, matrix):
    """Fraction of ALL spots with count > 0, per gene symbol (var index)."""
    import anndata as ad
    from scipy import sparse
    data = ad.read_h5ad(h5ad)
    X = data.layers[matrix[6:]] if matrix.startswith('layer:') else (data.raw.X if matrix == 'raw' else data.X)
    nnz = np.asarray((X > 0).sum(axis=0)).ravel() if sparse.issparse(X) else (np.asarray(X) > 0).sum(axis=0)
    names = data.var_names.astype(str)
    check(names.is_unique, 'Duplicate gene symbols in h5ad')
    return pd.Series(nnz / data.n_obs, index=names), data.n_obs


def sharing(pairs):
    """Counts for one set of pairs: size, distinct entities, and shared ligands/receptors."""
    if pairs.empty:
        return dict(pairs=0, ligands=0, receptors=0, pairs_sharing_ligand=0, pairs_sharing_receptor=0,
                    pairs_sharing_either=0, ligands_in_2plus=0, receptors_in_2plus=0)
    lc, rc = pairs.ligand.map(pairs.ligand.value_counts()), pairs.receptor.map(pairs.receptor.value_counts())
    return dict(pairs=len(pairs), ligands=pairs.ligand.nunique(), receptors=pairs.receptor.nunique(),
                pairs_sharing_ligand=int((lc > 1).sum()), pairs_sharing_receptor=int((rc > 1).sum()),
                pairs_sharing_either=int(((lc > 1) | (rc > 1)).sum()),
                ligands_in_2plus=int((pairs.ligand.value_counts() > 1).sum()),
                receptors_in_2plus=int((pairs.receptor.value_counts() > 1).sum()))


def annotate(it, detected, receptor_tf, tf_tg):
    """Add step flags (a, b, c) and target-set columns to every Secreted/Contact pair."""
    from evaluate_targets import target_sets
    p = it[it.annotation.isin(CLASSES)].copy()
    p['step_a'] = ~p.ligand_complex & ~p.receptor_complex
    p['ligand_detected'] = p.ligand.map(detected)
    p['receptor_detected'] = p.receptor.map(detected)
    p['ligand_in_data'] = p.ligand.isin(detected.index)
    p['receptor_in_data'] = p.receptor.isin(detected.index)
    p['step_b'] = p.step_a & (p.ligand_detected >= MIN_DETECTED) & (p.receptor_detected >= MIN_DETECTED)
    pool = set(p.loc[p.step_b, 'ligand']) | set(p.loc[p.step_b, 'receptor'])
    receptors = sorted(set(p.loc[p.step_b, 'receptor']))
    det = detected.to_dict()
    ident = {r: r for r in receptors}
    main = target_sets(receptors, ident, receptor_tf, tf_tg, det, pool, TOP_TF, MIN_TARGETS, MIN_DETECTED)
    main = main.set_index('receptor')
    p['receptor_tf_link'] = p.receptor.map(main.n_linked_tfs).fillna(0).astype(int) > 0
    p['n_targets'] = p.receptor.map(main.n_targets)
    p['target_status'] = p.receptor.map(main.status)
    p['step_c'] = p.step_b & (p.target_status == 'evaluated')
    # sensitivity: exclude only this pair's own ligand and receptor
    own = []
    for row in p.itertuples(index=False):
        if not row.step_b:
            own.append(np.nan)
            continue
        t = target_sets([row.receptor], {row.receptor: row.receptor}, receptor_tf, tf_tg, det,
                        {row.ligand, row.receptor}, TOP_TF, MIN_TARGETS, MIN_DETECTED).iloc[0]
        own.append(t.n_targets)
    p['n_targets_own_exclusion'] = own
    p['step_c_own_exclusion'] = p.step_b & p.receptor_tf_link & (p.n_targets_own_exclusion >= MIN_TARGETS)
    return p, main.reset_index()


def step_table(p):
    rows = []
    for cls in CLASSES + ['both classes']:
        q = p if cls == 'both classes' else p[p.annotation == cls]
        a = q[q.step_a]
        rows.append(dict(cls=cls, step='0 all', **sharing(q)))
        rows.append(dict(cls=cls, step='a single-gene L and R', **sharing(a)))
        b = q[q.step_b]
        rows.append(dict(cls=cls, step='b both detected >=5% spots', **sharing(b)))
        rows.append(dict(cls=cls, step='c1 receptor has TF link', **sharing(b[b.receptor_tf_link])))
        rows.append(dict(cls=cls, step='c >=5 targets (top-10 TF)', **sharing(q[q.step_c])))
        rows.append(dict(cls=cls, step='c sensitivity: own L/R excluded', **sharing(q[q.step_c_own_exclusion])))
    return pd.DataFrame(rows)


def complex_table(it):
    rows = []
    for cls in CLASSES:
        q = it[it.annotation == cls]
        rows.append(dict(cls=cls, single_gene=int((~q.ligand_complex & ~q.receptor_complex).sum()),
                         complex_ligand_only=int((q.ligand_complex & ~q.receptor_complex).sum()),
                         complex_receptor_only=int((~q.ligand_complex & q.receptor_complex).sum()),
                         complex_both=int((q.ligand_complex & q.receptor_complex).sum()),
                         total=len(q)))
    return pd.DataFrame(rows)


def conflicts(it, p, ccc):
    """CellChat single-gene pairs (all annotations) matched to CCCvelo by upper-case symbols."""
    single = it[~it.ligand_complex & ~it.receptor_complex].copy()
    single['L'], single['R'] = single.ligand.str.upper(), single.receptor.str.upper()
    table = ccc.rename(columns={'type': 'cccvelo_type', 'origin': 'cccvelo_origin'})
    same = single.merge(table.rename(columns={'ligand': 'L', 'receptor': 'R'}), on=['L', 'R'])
    # CCCvelo lists the pair with ligand and receptor swapped (e.g. NLGN1->NRXN1 vs Nrxn1->Nlgn1)
    swapped = single.merge(table.rename(columns={'ligand': 'R', 'receptor': 'L'}), on=['L', 'R'])
    swapped = swapped[swapped.L != swapped.R]  # homophilic pairs: swapped == same
    m = pd.concat([same.assign(orientation='same'), swapped.assign(orientation='swapped')], ignore_index=True)
    m['expected'] = m.annotation.map(CCC_TYPE)
    m['relation'] = np.where(m.expected.isna(), 'cellchat_class_without_cccvelo_counterpart',
                             np.where(m.expected == m.cccvelo_type, 'agree', 'conflict'))
    flags = p.set_index(['ligand', 'receptor'])[['step_b', 'step_c']]
    m = m.join(flags, on=['ligand', 'receptor'])
    cols = ['ligand', 'receptor', 'annotation', 'pathway_name', 'evidence', 'cccvelo_type', 'cccvelo_origin',
            'orientation', 'relation', 'step_b', 'step_c']
    return m[cols].sort_values(['relation', 'orientation', 'annotation', 'pathway_name', 'ligand', 'receptor'])


def report(args):
    from evaluate_targets import load_db
    from stage1 import digest, load_lr, write_json
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    it, meta = load_cellchat(args.cellchat)
    receptor_tf, tf_tg, tmeta = load_db(args.scseqcomm)
    detected, n_spots = detected_fractions(args.h5ad, args.matrix)
    p, receptors = annotate(it, detected, receptor_tf, tf_tg)
    steps, cx = step_table(p), complex_table(it)
    ccc = load_lr(args.lr_dir)
    conf = conflicts(it, p, ccc)
    panel = pd.read_csv(args.panel, dtype=str) if args.panel else None
    out.mkdir(parents=True)
    steps.to_csv(out/'steps.csv', index=False)
    cx.to_csv(out/'complexes.csv', index=False)
    p.drop(columns=['ligand_complex', 'receptor_complex']).to_csv(out/'pairs.csv', index=False)
    receptors.to_csv(out/'receptor_targets.csv', index=False)
    conf.to_csv(out/'cellchat_vs_cccvelo.csv', index=False)
    summary = dict(status='complete', h5ad=str(Path(args.h5ad).resolve()), h5ad_sha256=digest(args.h5ad),
        matrix=args.matrix, n_spots=int(n_spots), detection='count > 0 in >= 5% of all spots',
        cellchat=dict((k, meta[k]) for k in ['source', 'commit', 'file_last_changed_commit', 'raw_sha256',
                                             'annotation_counts']),
        scseqcomm=tmeta['source'], cccvelo_lr=f'CCCvelo pinned tables ({args.lr_dir})',
        rules=dict(top_tf=TOP_TF, min_targets=MIN_TARGETS, min_detected=MIN_DETECTED,
                   target_exclusion='all single-gene L/R genes of step-b pairs (both classes); '
                                    'sensitivity: own L/R only'),
        relation_counts=conf.relation.value_counts().to_dict(),
        note='Inventory only; no panel selected.')
    if panel is not None:
        cc = it[['ligand', 'receptor', 'annotation', 'pathway_name']].assign(
            L=it.ligand.str.upper(), R=it.receptor.str.upper())
        rows = []
        for row in panel.itertuples(index=False):
            same = cc[(cc.L == row.ligand) & (cc.R == row.receptor)]
            swap = cc[(cc.L == row.receptor) & (cc.R == row.ligand)]
            hit, orient = (same, 'same') if len(same) else (swap, 'swapped') if len(swap) else (None, 'absent')
            rows.append(dict(ligand=row.ligand, receptor=row.receptor, cccvelo_type=row.type,
                             cellchat_orientation=orient,
                             cellchat_annotation=None if hit is None else hit.annotation.iloc[0],
                             cellchat_pathway=None if hit is None else hit.pathway_name.iloc[0]))
        pd.DataFrame(rows).to_csv(out/'current_panel_in_cellchat.csv', index=False)
    write_json(out/'inventory.json', summary)
    pd.set_option('display.width', 220)
    print(cx.to_string(index=False))
    print(steps.to_string(index=False))
    print(conf.relation.value_counts().to_string())


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fetch'); f.add_argument('--output', default='data/cellchatdb')
    r = sub.add_parser('report')
    r.add_argument('--h5ad', default='data/mousebrain_bin60.h5ad')
    r.add_argument('--matrix', default='layer:count')
    r.add_argument('--cellchat', default='data/cellchatdb')
    r.add_argument('--scseqcomm', default='data/scseqcomm')
    r.add_argument('--lr-dir', default='data/cccvelo_lr')
    r.add_argument('--panel', help='current lr.csv, to show its CellChat annotation')
    r.add_argument('--output', required=True)
    args = p.parse_args()
    fetch(args.output) if args.command == 'fetch' else report(args)


if __name__ == '__main__':
    main()
