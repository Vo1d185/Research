"""Target-gene evaluation of inferred communication (the COMMOT-style check).

Idea (COMMOT, Nat Methods 2023, Methods "Evaluation metrics"): spots that receive
more signal through a receptor should express more of its known downstream genes.
COMMOT used scSeqComm target libraries (mouse: TF_TG_TRRUSTv2_RegNetwork_High_mouse)
and Spearman correlation between received signal and mean target expression.

  fetch-db  download the pinned scSeqComm mouse tables (GPL-3, not redistributed)
            and convert them to receptor_tf.csv / tf_tg.csv
  run       score any number of result directories (research.py, stage1 runs,
            commot_baseline.py) against the same target sets

Pre-specified rules (change only with a documented reason, never after seeing results):
  * receptor -> TF: scSeqComm KEGG+REACTOME mouse PPR scores; a TF's weight is its
    maximum PPR over pathways; keep the top-k TFs (k=10 main, 5 and 20 sensitivity).
    COMMOT's exact LR->target rule is in its Zenodo reproduction code (10.5281/
    zenodo.7272562), not verified here; align when available.
  * targets: union of those TFs' targets, minus every LR-panel gene, kept if
    detected in >= 5% of ROI spots; a receptor needs >= 5 targets.
  * target score per spot: mean normalized expression (COMMOT); z-scored mean as sensitivity.
  * received signal: sum of T over sender spots, per LR pair and per receptor.
  * metrics: Spearman(received, target); receptor-only Spearman (receptor expression
    instead of received signal); partial Spearman controlling receptor expression
    (+ optional cell-type covariate). No p-values: spots are spatially autocorrelated.

Spearman is rank-based: if PPI scores only rescale an LR pair's whole plan, its
per-pair Spearman cannot change. rank_corr_vs_no_ppi reports this directly.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from research import check

SCSEQCOMM_COMMIT = '95c9bb251cd94eb8d38db41086e71c1bf46f625a'
SCSEQCOMM_RAW = f'https://gitlab.com/sysbiobig/scseqcomm/-/raw/{SCSEQCOMM_COMMIT}/data/'
SCSEQCOMM_FILES = {
    'TF_TG_TRRUSTv2_RegNetwork_High_mouse.rda': '7b8299c30ce4908d8f90f345afd4ed9360ac0f2a9feb186150682f0444ed36e3',
    'TF_PPR_KEGG_mouse.rda': '73a4802219326a2b352e2a4c35c3519b828eb8d9d356327c68b97abfa37ff622',
    'TF_PPR_REACTOME_mouse.rda': '0d5c826045e49b1386fd6e684b20af9e37c7707ecd70ce433415bfcb6c0259af',
}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False, default=str))


# ------------------------------------------------------------------ database
def fetch_db(output):
    import rdata
    from stage1 import download
    out = Path(output)
    check(not (out/'db.json').exists(), f'{out} already has a converted database')
    import time
    from urllib.error import HTTPError
    for name, sha in SCSEQCOMM_FILES.items():
        for attempt in range(4):  # GitLab occasionally rate-limits anonymous raw downloads
            try:
                download(SCSEQCOMM_RAW + name, out/'raw'/name, expected=sha)
                break
            except HTTPError as error:
                if error.code not in (429, 500, 502, 503, 504) or attempt == 3:
                    raise
                time.sleep(5 * 2**attempt)
    tf_tg = rdata.read_rda(out/'raw'/'TF_TG_TRRUSTv2_RegNetwork_High_mouse.rda')
    tf_tg = next(iter(tf_tg.values()))
    edges = pd.DataFrame([(str(tf), str(t)) for tf, targets in tf_tg.items() for t in np.atleast_1d(targets)],
                         columns=['tf', 'target']).drop_duplicates()
    frames = []
    for name, source in [('TF_PPR_KEGG_mouse.rda', 'KEGG'), ('TF_PPR_REACTOME_mouse.rda', 'REACTOME')]:
        frame = next(iter(rdata.read_rda(out/'raw'/name).values()))
        frames.append(pd.DataFrame(dict(receptor=frame.receptor.astype(str), pathway=frame.pathway.astype(str),
                                        tf=frame.tf.astype(str), weight=frame.tf_PPR.astype(float), source=source)))
    receptor_tf = pd.concat(frames, ignore_index=True)
    edges.to_csv(out/'tf_tg.csv', index=False)
    receptor_tf.to_csv(out/'receptor_tf.csv', index=False)
    write_json(out/'db.json', dict(status='complete', source=f'scSeqComm {SCSEQCOMM_COMMIT}',
        license='GPL-3 (downloaded, not redistributed)', raw_sha256=SCSEQCOMM_FILES,
        weight='personalized PageRank of TF seeded at the receptor (scSeqComm TF_PPR)',
        tf_tg_edges=len(edges), tfs=int(edges.tf.nunique()), receptor_tf_rows=len(receptor_tf),
        receptors=int(receptor_tf.receptor.nunique()),
        converted_sha256={f: digest(out/f) for f in ['tf_tg.csv', 'receptor_tf.csv']}))
    print(f'Converted scSeqComm mouse tables -> {out}')


def load_db(folder):
    folder = Path(folder)
    meta = json.loads((folder/'db.json').read_text())
    check(meta.get('status') == 'complete', 'Incomplete target database')
    for name, sha in meta['converted_sha256'].items():
        check(digest(folder/name) == sha, f'{name} changed after conversion')
    receptor_tf = pd.read_csv(folder/'receptor_tf.csv', dtype={'receptor': str, 'tf': str})
    tf_tg = pd.read_csv(folder/'tf_tg.csv', dtype=str)
    check({'receptor', 'tf'} <= set(receptor_tf) and {'tf', 'target'} <= set(tf_tg),
          'receptor_tf.csv needs receptor,tf[,weight]; tf_tg.csv needs tf,target')
    return receptor_tf, tf_tg, meta


def target_sets(receptors, symbol, receptor_tf, tf_tg, detected, panel_genes, top_tf,
                min_targets=5, min_detected=0.05):
    """receptors: LR symbols; symbol: LR symbol -> data symbol; detected: data symbol -> fraction."""
    targets_of = tf_tg.groupby('tf').target.apply(set).to_dict()
    weighted = 'weight' in receptor_tf and receptor_tf.weight.notna().any()
    rows = []
    for receptor in receptors:
        data_symbol = symbol[receptor]
        links = receptor_tf[receptor_tf.receptor == data_symbol]
        if weighted:
            ranked = links.groupby('tf').weight.max().reset_index().sort_values(['weight', 'tf'],
                                                                                 ascending=[False, True])
            tfs = ranked.tf.head(top_tf).tolist()
        else:  # unweighted database: no defensible top-k, use all linked TFs
            tfs = sorted(links.tf.unique())
        candidates = set().union(*[targets_of.get(tf, set()) for tf in tfs]) if tfs else set()
        kept = sorted(g for g in candidates - panel_genes if detected.get(g, 0) >= min_detected)
        status = ('no_receptor_tf_link' if not len(links) else
                  'too_few_targets' if len(kept) < min_targets else 'evaluated')
        rows.append(dict(receptor=receptor, data_symbol=data_symbol, top_tf=top_tf if weighted else 'all',
                         n_linked_tfs=int(links.tf.nunique()), tfs=';'.join(tfs),
                         n_candidate_targets=len(candidates), n_targets=len(kept), status=status,
                         targets=';'.join(kept)))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------ statistics
def spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(np.corrcoef(rankdata(x), rankdata(y))[0, 1])


def partial_spearman(x, y, controls):
    """Rank-based partial correlation; controls: list of 1-D arrays (ranked) or 2-D dummies."""
    x, y = rankdata(x), rankdata(y)
    Z = [np.ones(len(x))]
    for c in controls:
        c = np.asarray(c, float)
        Z.extend([rankdata(c)] if c.ndim == 1 else list(c.T))
    Z = np.column_stack(Z)
    rx = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    ry = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    if np.ptp(rx) < 1e-12 or np.ptp(ry) < 1e-12:
        return np.nan
    return float(np.corrcoef(rx, ry)[0, 1])


def target_scores(X, genes, targets, kind):
    index = {g: i for i, g in enumerate(genes)}
    sub = X[:, [index[g] for g in targets]].astype(float)
    if kind == 'mean':
        return sub.mean(axis=1)
    sd = sub.std(axis=0)
    keep = sd > 0
    return ((sub[:, keep] - sub[:, keep].mean(axis=0)) / sd[keep]).mean(axis=1)


# ------------------------------------------------------------------ evaluation
def received_signal(T, ligands, receptors, n, pairs):
    li, ri = {v: i for i, v in enumerate(ligands)}, {v: i for i, v in enumerate(receptors)}
    per_pair = {}
    for row in pairs.itertuples(index=False):
        i, j = li[row.ligand], ri[row.receptor]
        per_pair[(row.ligand, row.receptor)] = T[i*n:(i+1)*n, j*n:(j+1)*n].sum(axis=0)
    per_receptor = {}
    for (ligand, receptor), signal in per_pair.items():
        per_receptor[receptor] = per_receptor.get(receptor, 0) + signal
    return per_pair, per_receptor


def load_results(folder, input_dir):
    folder = Path(folder)
    manifest = json.loads((folder/'manifest.json').read_text())
    check(manifest.get('status') == 'complete', f'{folder}: incomplete results')
    for name in ['expression.csv', 'spatial.csv', 'lr.csv']:
        check(manifest['input_sha256'].get(name) == digest(Path(input_dir)/name),
              f'{folder}: {name} differs from --input; results were computed on other inputs')
    summary = pd.read_csv(folder/'summary.csv')
    pairs = pd.read_csv(folder/'used_lr.csv', dtype=str)
    return manifest, summary, pairs


def evaluate(args):
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    inp, genes_dir = Path(args.input), Path(args.genes)
    export = json.loads((genes_dir/'gene_export.json').read_text())
    check(export.get('status') == 'complete' and digest(genes_dir/'gene_expression.npz') ==
          export['gene_expression_sha256'], 'Gene export incomplete or changed')
    with np.load(genes_dir/'gene_expression.npz', allow_pickle=False) as z:
        X, spots, genes = z['X'], z['spots'].astype(str), z['genes'].astype(str)
    spatial = pd.read_csv(inp/'spatial.csv', index_col=0)
    check(list(spots) == list(spatial.index.astype(str)), 'Gene export spots differ from --input')
    gene_table = pd.read_csv(genes_dir/'genes.csv', dtype={'data_symbol': str, 'lr_symbol': str})
    symbol = dict(zip(gene_table.lr_symbol.dropna(), gene_table.data_symbol[gene_table.lr_symbol.notna()]))
    detected = dict(zip(gene_table.data_symbol, gene_table.detected_fraction))
    expression = pd.read_csv(inp/'expression.csv', index_col=0)
    expression.index = expression.index.astype(str)
    expression = expression.loc[spots]
    covariate = []
    if args.covariate:
        meta = pd.read_csv(inp/'cell_meta.csv', index_col=0)
        meta.index = meta.index.astype(str)
        check(args.covariate in meta, f'{args.covariate} not in cell_meta.csv')
        dummies = pd.get_dummies(meta.loc[spots, args.covariate].astype(str), drop_first=True)
        covariate = [dummies.to_numpy(float)] if dummies.shape[1] else []
    receptor_tf, tf_tg, db_meta = load_db(args.db)
    results = [(Path(r),) + load_results(r, inp) for r in args.results]
    panel = pd.read_csv(inp/'lr.csv', dtype=str)
    panel_genes = {symbol[g] for g in set(panel.ligand) | set(panel.receptor) if g in symbol}
    receptors = sorted(set(panel.receptor))
    missing = [r for r in receptors if r not in symbol]
    check(not missing, f'Receptors without data symbol: {missing}')

    sets = pd.concat([target_sets(receptors, symbol, receptor_tf, tf_tg, detected, panel_genes, k,
                                  args.min_targets, args.min_detected) for k in args.top_tf])
    sets = sets.drop_duplicates(['receptor', 'top_tf'])
    primary_k = sets.top_tf.iloc[0]
    scores = {}
    for row in sets[sets.status == 'evaluated'].itertuples(index=False):
        for kind in args.target_score:
            scores[(row.receptor, row.top_tf, kind)] = target_scores(X, genes, row.targets.split(';'), kind)
    receptor_expr = {r: expression[r].to_numpy(float) for r in receptors}

    rows = []
    for folder, manifest, summary, pairs in results:
        refs = {}
        for case in summary.itertuples(index=False):
            with np.load(folder/f'{case.case}.npz', allow_pickle=False) as z:
                T = z['T']
                check(list(z['spots'].astype(str)) == list(spots), f'{folder}/{case.case}: spot order differs')
                ligands, recs = list(z['ligands'].astype(str)), list(z['receptors'].astype(str))
            per_pair, per_receptor = received_signal(T, ligands, recs, len(spots), pairs)
            key = (str(folder), case.family, None if pd.isna(case.r0) else float(case.r0))
            units = [('pair', f'{l}->{r}', r, s, sum(pairs.receptor == r)) for (l, r), s in per_pair.items()]
            units += [('receptor', r, r, s, sum(pairs.receptor == r)) for r, s in per_receptor.items()]
            if case.score_mode == 'no_ppi':
                refs[key] = {(u[0], u[1]): u[3] for u in units}
            for level, unit, receptor, signal, n_ligands in units:
                reference = refs.get(key, {}).get((level, unit))
                for (r, k, kind), target in scores.items():
                    if r != receptor:
                        continue
                    rows.append(dict(results=str(folder), method=manifest.get('method', 'research.py'),
                        case=case.case, family=case.family, r0=case.r0,
                        r0_h=getattr(case, 'dis_thr_h', None) if case.family == 'commot_euc'
                        else (None if pd.isna(case.r0) else case.r0/manifest['h']),
                        score_mode=case.score_mode, permuted_type=getattr(case, 'permuted_type', None),
                        level=level, unit=unit, receptor=receptor, n_ligands_of_receptor=int(n_ligands),
                        top_tf=k, target_score=kind, total_received=float(signal.sum()),
                        spearman=spearman(signal, target),
                        partial_spearman=partial_spearman(signal, target, [receptor_expr[receptor]] + covariate),
                        receptor_only_spearman=spearman(receptor_expr[receptor], target),
                        rank_corr_vs_no_ppi=spearman(signal, reference) if reference is not None else np.nan))
    metrics = pd.DataFrame(rows)
    check(len(metrics) > 0, 'No receptor has an evaluable target set; see target_sets.csv')
    out.mkdir(parents=True)
    sets.to_csv(out/'target_sets.csv', index=False)
    metrics.to_csv(out/'case_metrics.csv', index=False)
    main = metrics[(metrics.top_tf.astype(str) == str(primary_k)) & (metrics.target_score == args.target_score[0])]
    method = (main[main.score_mode == 'no_ppi'].groupby(['results', 'method', 'family', 'r0', 'r0_h', 'level'],
              dropna=False)[['spearman', 'partial_spearman', 'receptor_only_spearman']].mean().reset_index())
    method.to_csv(out/'method_summary.csv', index=False)
    comparison, permutation = ppi_tables(main)
    comparison.to_csv(out/'ppi_comparison.csv', index=False)
    permutation.to_csv(out/'permutation_summary.csv', index=False)
    plot(method, main, out/'target_eval.png')
    write_json(out/'eval_manifest.json', dict(status='complete', input=str(inp.resolve()),
        gene_export_sha256=export['gene_expression_sha256'], database=db_meta,
        results=[str(r[0]) for r in results], primary_top_tf=primary_k, top_tf=args.top_tf,
        target_score=args.target_score, min_targets=args.min_targets, min_detected=args.min_detected,
        covariate=args.covariate, receptors_evaluated=sorted(set(metrics.receptor)),
        receptors_not_evaluated=sets[sets.status != 'evaluated'][['receptor', 'top_tf', 'status']].to_dict('records'),
        notes=['Effect sizes only; spatial autocorrelation invalidates naive p-values.',
               'rank_corr_vs_no_ppi == 1 means PPI only rescaled that unit: Spearman cannot change.',
               'Receptor-level units sum all ligands of a receptor; they can change when scores differ.']))
    print(f'Evaluated {metrics.receptor.nunique()} receptors -> {out}')


def ppi_tables(main):
    """Real scores vs S=1, the type-mean constant, and within-type permutations."""
    comp, perm = [], []
    keys = ['results', 'family', 'r0', 'level', 'unit']
    for key, group in main.groupby(keys, dropna=False):
        base = group[group.score_mode == 'no_ppi']
        real = group[group.score_mode == 'ppi']
        if len(base) != 1 or len(real) != 1:
            continue
        base, real = base.iloc[0], real.iloc[0]
        for row in group[group.score_mode != 'no_ppi'].itertuples(index=False):
            comp.append(dict(zip(keys, key), score_mode=row.score_mode, permuted_type=row.permuted_type,
                             delta_spearman=row.spearman - base.spearman,
                             delta_partial_spearman=row.partial_spearman - base.partial_spearman,
                             rank_corr_vs_no_ppi=row.rank_corr_vs_no_ppi))
        shuffled = group[group.score_mode.str.startswith(('perm_', 'shuffled_'))]
        for ptype, sub in shuffled.groupby(shuffled.permuted_type.fillna('both')):
            for metric in ['spearman', 'partial_spearman']:
                values = sub[metric].dropna().to_numpy()
                if not len(values) or np.isnan(real[metric]):
                    continue
                perm.append(dict(zip(keys, key), permuted_type=ptype, metric=metric, real=real[metric],
                    n_permutations=len(values), permutation_mean=float(values.mean()),
                    # (1 + #perm >= real)/(1 + #perm): exact enumeration, small panels => coarse
                    p_ge=float((1 + (values >= real[metric] - 1e-12).sum()) / (1 + len(values)))))
    return pd.DataFrame(comp), pd.DataFrame(perm)


def plot(method, main, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from analyze import PALETTE, style_axis
    pair = method[method.level == 'pair'].copy()
    pair['label'] = [f'{f} {r:g}h' if pd.notna(r) else f for f, r in zip(pair.family, pair.r0_h)]
    # The same S=1 design can appear in several result dirs (e.g. S=1 run and PPI run):
    # identical metrics are shown once; differing ones are labelled with their result dir.
    pair = pair.drop_duplicates(['label', 'spearman', 'partial_spearman'])
    repeated = pair.label.duplicated(keep=False)
    pair.loc[repeated, 'label'] += ' [' + pair.loc[repeated, 'results'].map(lambda p: '/'.join(Path(p).parts[-3:])) + ']'
    pair = pair.sort_values('partial_spearman')
    has_ppi = (main.score_mode == 'ppi').any()
    fig, axes = plt.subplots(1, 2 if has_ppi else 1, figsize=(12 if has_ppi else 6.5, 0.45*len(pair) + 2),
                             constrained_layout=True, squeeze=False)
    ax = axes[0, 0]
    y = np.arange(len(pair))
    ax.scatter(pair.spearman, y, color=PALETTE['series2'], marker='o', label='Spearman', zorder=3, s=30)
    ax.scatter(pair.partial_spearman, y, color=PALETTE['series1'], marker='D', label='partial (| receptor)', zorder=3, s=30)
    receptor_only = pair.receptor_only_spearman.mean()
    ax.axvline(receptor_only, color=PALETTE['muted'], ls='--', lw=1, label='receptor expression only')
    ax.set_yticks(y, pair.label)
    ax.set_xlabel('Mean correlation with target genes (S = 1, per LR pair)')
    ax.set_title('Distance designs and COMMOT baseline')
    ax.legend(fontsize=8, frameon=False, loc='lower right')
    style_axis(ax)
    if has_ppi:
        ax = axes[0, 1]
        sub = main[(main.level == 'pair') & main.results.isin(main[main.score_mode == 'ppi'].results)]
        units = sorted(sub.unit.unique())
        for k, unit in enumerate(units):
            u = sub[sub.unit == unit]
            perm = u[u.score_mode.str.startswith(('perm_', 'shuffled_'))]
            ax.scatter(perm.partial_spearman, np.full(len(perm), k), color=PALETTE['muted'], s=10,
                       alpha=.6, label='permuted scores' if k == 0 else None, zorder=2)
            for mode, color, marker in [('no_ppi', PALETTE['series2'], 'o'), ('ppi_constant', PALETTE['series3'], 's'),
                                        ('ppi', PALETTE['series1'], 'D')]:
                m = u[u.score_mode == mode]
                ax.scatter(m.partial_spearman, np.full(len(m), k), color=color, marker=marker, s=34, zorder=3,
                           label={'no_ppi': 'S = 1', 'ppi_constant': 'type-mean constant',
                                  'ppi': 'real PLM scores'}[mode] if k == 0 else None)
        ax.set_yticks(np.arange(len(units)), units)
        ax.set_xlabel('Partial Spearman with target genes')
        ax.set_title('Real scores vs controls (same distance design)')
        ax.legend(fontsize=8, frameon=False, loc='lower right')
        style_axis(ax)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fetch-db'); f.add_argument('--output', default='data/scseqcomm')
    r = sub.add_parser('run')
    r.add_argument('--input', required=True, help='the S=1 prepared input (stage1_input_01)')
    r.add_argument('--genes', required=True, help='output of stage1.py export-genes')
    r.add_argument('--db', default='data/scseqcomm')
    r.add_argument('--results', nargs='+', required=True, help='result dirs (with summary.csv, case_*.npz)')
    r.add_argument('--output', required=True)
    r.add_argument('--top-tf', type=int, nargs='+', default=[10, 5, 20], help='first value is the main analysis')
    r.add_argument('--target-score', nargs='+', choices=['mean', 'zmean'], default=['mean', 'zmean'])
    r.add_argument('--min-targets', type=int, default=5)
    r.add_argument('--min-detected', type=float, default=0.05)
    r.add_argument('--covariate', help='column of cell_meta.csv, e.g. a cell-type annotation')
    return p


if __name__ == '__main__':
    arguments = parser().parse_args()
    fetch_db(arguments.output) if arguments.command == 'fetch-db' else evaluate(arguments)
