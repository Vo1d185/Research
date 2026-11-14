"""Pre-registered, single, time-boxed check of PLM-interact scores (commot env).

Criteria are fixed in diagnostics/README.md ("预先登记的最后一次尝试") and were
committed BEFORE any Leakage-Free score or any double-centered score was computed.
Do not change thresholds, transforms, candidates or inputs after seeing results.

Usage: one --checkpoint per model, each with its score files for
  (1) the human positive-control run, (2) the mouse panel full length,
  (3) the mouse panel extracellular region:
    --checkpoint NAME CONTROLS_SCORES PANEL_FULL_SCORES PANEL_ECD_SCORES
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from proteins import auroc, digest, write_json  # noqa: E402
from research import check  # noqa: E402

CLIP = 1e-6
AUROC_MIN = 0.85
FIRST_MIN = 10
CANDIDATES = {('LeakageFree', 'raw'), ('V11', 'double_centered'),
              ('V12', 'double_centered'), ('LeakageFree', 'double_centered')}


def pair_table(jobs_path, scores_path):
    """Both input orders per pair -> ligand, receptor, set, p_LR, p_RL."""
    jobs = pd.read_csv(jobs_path, usecols=lambda c: c in {'pair_id', 'set', 'ligand', 'receptor', 'order'})
    scores = pd.read_csv(scores_path)
    check(set(scores.pair_id) == set(jobs.pair_id), f'{scores_path}: pair IDs differ from {jobs_path}')
    m = jobs.merge(scores[['pair_id', 'score']], on='pair_id', validate='one_to_one')
    lr = m[m.order == 'LR'][['ligand', 'receptor', 'set', 'score']].rename(columns={'score': 'p_LR'})
    rl = m[m.order == 'RL'][['ligand', 'receptor', 'score']].rename(columns={'score': 'p_RL'})
    t = lr.merge(rl, on=['ligand', 'receptor'], validate='one_to_one')
    check(len(t) * 2 == len(m), f'{scores_path}: every pair needs both orders')
    return t


def logit(p):
    p = np.clip(np.asarray(p, float), CLIP, 1 - CLIP)
    return np.log(p / (1 - p))


def transform(t, kind):
    """raw: mean probability of both orders. double_centered: mean logit of both
    orders, then z_ij = x_ij - row_mean_i - col_mean_j + grand_mean over the FULL
    ligand x receptor grid of the set (all pairs, cognate included)."""
    t = t.copy()
    if kind == 'raw':
        t['value'] = (t.p_LR + t.p_RL) / 2
        return t
    check(kind == 'double_centered', f'Unknown transform {kind}')
    t['x'] = (logit(t.p_LR) + logit(t.p_RL)) / 2
    grid = t.pivot(index='ligand', columns='receptor', values='x')
    check(not grid.isna().any().any(), 'Double centering needs the complete ligand x receptor grid')
    z = grid.sub(grid.mean(axis=1), axis=0).sub(grid.mean(axis=0), axis=1) + grid.to_numpy().mean()
    t['value'] = [z.loc[l, r] for l, r in zip(t.ligand, t.receptor)]
    return t


def controls_metrics(t, known):
    """AUROC (cognate vs unannotated cross, documented cross-reactive excluded) and
    number of ligands whose cognate receptor scores strictly above every competitor
    (competitors exclude documented cross-reactive receptors of that ligand)."""
    is_known = pd.Series([(l, r) in known for l, r in zip(t.ligand, t.receptor)], index=t.index)
    cog = t[t.set == 'cognate']
    cross = t[(t.set == 'cross') & ~is_known]
    ranks = []
    for row in cog.itertuples(index=False):
        comp = t[(t.ligand == row.ligand) & (t.receptor != row.receptor) & ~is_known]
        check(len(comp) > 0, f'{row.ligand}: no competing receptors')
        ranks.append(dict(ligand=row.ligand, receptor=row.receptor, value=row.value,
                          best_competitor=comp.loc[comp.value.idxmax(), 'receptor'],
                          best_competitor_value=comp.value.max(),
                          rank=int((comp.value >= row.value).sum()) + 1,
                          first=bool(row.value > comp.value.max())))
    return auroc(cog.value, cross.value), pd.DataFrame(ranks)


def direction(t, kind):
    """Per cognate pair: +1/-1/0. double_centered: sign(z). raw: sign(S - median S
    over all pairs of the panel grid in the same region)."""
    ref = 0.0 if kind == 'double_centered' else float(t.value.median())
    cog = t[t.set == 'cognate']
    return pd.Series(np.sign(cog.value.to_numpy() - ref), index=list(zip(cog.ligand, cog.receptor)))


def evaluate(checkpoints, controls_jobs, full_jobs, ecd_jobs, known):
    rows, ranks, dirs = [], [], []
    for name, (cpath, fpath, epath) in checkpoints.items():
        tables = {'controls': pair_table(controls_jobs, cpath), 'full': pair_table(full_jobs, fpath),
                  'ecd': pair_table(ecd_jobs, epath)}
        for set_name, frame in [('full', tables['full']), ('ecd', tables['ecd'])]:
            check(set(frame.set) <= {'cognate', 'panel_cross'}, f'{set_name}: unexpected pair sets')
        for kind in ['raw', 'double_centered']:
            c = transform(tables['controls'], kind)
            auc, r = controls_metrics(c, known)
            d_full = direction(transform(tables['full'], kind), kind)
            d_ecd = direction(transform(tables['ecd'], kind), kind)
            check(list(d_full.index) == list(d_ecd.index), 'Full and ECD cognate pairs differ')
            consistent = (d_full == d_ecd) & (d_full != 0)
            n_first = int(r['first'].sum())
            rows.append(dict(checkpoint=name, transform=kind, candidate=(name, kind) in CANDIDATES,
                auroc=auc, n_first=n_first, n_ligands=len(r),
                direction_consistent=int(consistent.sum()), n_panel_cognate=len(consistent),
                pass_auroc=bool(auc >= AUROC_MIN), pass_rank=bool(n_first >= FIRST_MIN),
                pass_direction=bool(consistent.all())))
            rows[-1]['pass_all'] = rows[-1]['pass_auroc'] and rows[-1]['pass_rank'] and rows[-1]['pass_direction']
            ranks.append(r.assign(checkpoint=name, transform=kind))
            dirs.append(pd.DataFrame(dict(checkpoint=name, transform=kind,
                ligand=[k[0] for k in d_full.index], receptor=[k[1] for k in d_full.index],
                direction_full=d_full.to_numpy(), direction_ecd=d_ecd.to_numpy(),
                consistent=consistent.to_numpy())))
    return pd.DataFrame(rows), pd.concat(ranks, ignore_index=True), pd.concat(dirs, ignore_index=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--controls-jobs', required=True)
    p.add_argument('--panel-full-jobs', required=True)
    p.add_argument('--panel-ecd-jobs', required=True)
    p.add_argument('--known-cross', required=True)
    p.add_argument('--checkpoint', nargs=4, action='append', required=True,
                   metavar=('NAME', 'CONTROLS', 'FULL', 'ECD'))
    p.add_argument('--output', required=True)
    args = p.parse_args()
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    k = pd.read_csv(args.known_cross, dtype=str)
    known = set(zip(k.ligand, k.receptor))
    checkpoints = {c[0]: c[1:] for c in args.checkpoint}
    check({n for n, _ in CANDIDATES} <= set(checkpoints), 'Every candidate checkpoint must be supplied')
    summary, ranks, dirs = evaluate(checkpoints, args.controls_jobs, args.panel_full_jobs,
                                    args.panel_ecd_jobs, known)
    out.mkdir(parents=True)
    summary.to_csv(out/'criteria.csv', index=False)
    ranks.to_csv(out/'ligand_ranks.csv', index=False)
    dirs.to_csv(out/'panel_directions.csv', index=False)
    passed = summary[summary.candidate & summary.pass_all]
    write_json(out/'verdict.json', dict(status='complete',
        verdict='pass' if len(passed) else 'fail_all_candidates',
        passing_candidates=[f'{r.checkpoint}/{r.transform}' for r in passed.itertuples()],
        thresholds=dict(auroc_min=AUROC_MIN, first_min=FIRST_MIN, direction='all 8 cognate pairs', clip=CLIP),
        inputs_sha256={f'{n}/{role}': digest(path) for n, paths in checkpoints.items()
                       for role, path in zip(['controls', 'full', 'ecd'], paths)},
        note='Non-candidate rows (V11/V12 raw) are reference only and do not affect the verdict.'))
    pd.set_option('display.width', 200)
    print(summary.to_string(index=False, float_format=lambda v: f'{v:.3f}'))
    print('verdict:', 'pass ' + str(list(passed.checkpoint + '/' + passed.transform)) if len(passed)
          else 'all candidates fail')


if __name__ == '__main__':
    main()
