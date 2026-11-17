"""COMMOT-style target-gene evaluation values from existing results only (no new OT runs, no new evaluation).

Reads the completed stage-2 target evaluations (target_eval_stage2_{post,ant}/case_metrics.csv, written by
evaluate_targets.py) and summarises, per slide and condition, the per-receptor (and per-LR-pair) values of
  * raw Spearman (received signal vs mean target expression; receptor expression not removed),
  * partial Spearman (same, controlling receptor expression),
  * receptor-only Spearman (receptor expression vs mean target expression).
Conditions: C (COMMOT as published, tutorial cutoff), C726 (COMMOT, paper cutoff; supplementary), B, A at
r0 = 50/100/200/400 um, F (flat kernel), I (inverse kernel). Rules: main target rule (top-10 TFs) and the
COMMOT target rule (plan 8.8) where evaluated.

Statistics: median and mean over units; 95% percentile bootstrap intervals for both, resampling units with
replacement (1000 resamples, numpy.random.default_rng(1)), as in r0_validation.py.
Units: receptor level = receptors defined (non-NaN) in every listed condition of that slide and rule
("common set"); LR-pair level is reported with NaN kept out and, as in COMMOT's notebook, with NaN set to 0.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import stage2
from research import check, launch_info

REPLICATES, SEED = 1000, 1
ORDER = [('C', 0.0), ('C726', 0.0), ('B', 0.0), ('A', 50.0), ('A', 100.0), ('A', 200.0), ('A', 400.0),
         ('F', 0.0), ('I', 0.0)]
COMMON = [k for k in ORDER if k[0] != 'C726']     # the conditions the user listed; C726 is supplementary


def boot(values, seed=SEED):
    x = np.asarray(values, float)
    if len(x) == 0:
        return dict(n=0)
    idx = np.random.default_rng(seed).integers(0, len(x), size=(REPLICATES, len(x)))
    means, medians = x[idx].mean(axis=1), np.median(x[idx], axis=1)
    return dict(n=len(x), median=float(np.median(x)), median_ci_low=float(np.quantile(medians, .025)),
                median_ci_high=float(np.quantile(medians, .975)), mean=float(x.mean()),
                mean_ci_low=float(np.quantile(means, .025)), mean_ci_high=float(np.quantile(means, .975)),
                fraction_positive=float((x > 0).mean()))


def label(condition, r0):
    return f'A(r0={r0:g})' if condition == 'A' else condition


def load(ev):
    meta = json.loads((ev/'eval_manifest.json').read_text())
    check(meta.get('status') == 'complete', f'{ev}: incomplete evaluation')
    cases = pd.concat([stage2.condition_cases(r) for r in meta['results']], ignore_index=True)
    cases = cases[cases.condition.isin({c for c, _ in ORDER})].copy()
    cases['r0_key'] = cases.r0_um.fillna(0.0)
    check(not cases.duplicated(['condition', 'r0_key']).any(), f'{ev}: condition appears twice')
    m = pd.read_csv(ev/'case_metrics.csv', dtype={'top_tf': str})
    m = m[(m.target_score == meta['target_score'][0]) & (m.score_mode == 'no_ppi')]
    m = m.merge(cases[['results', 'case', 'condition', 'r0_key']], on=['results', 'case'], how='inner')
    lr = pd.read_csv(Path(meta['input'])/'lr.csv', dtype=str)
    return meta, m, lr


def receptor_tables(m, key, slide, rule):
    r = m[(m.level == 'receptor') & (m.top_tf == key)]
    check(not r.duplicated(['receptor', 'condition', 'r0_key']).any(), 'duplicate receptor rows')
    present = [k for k in ORDER if ((r.condition == k[0]) & (r.r0_key == k[1])).any()]
    wide = {v: r.set_index(['receptor', 'condition', 'r0_key'])[v].unstack(['condition', 'r0_key'])
            for v in ['spearman', 'partial_spearman']}
    listed = [k for k in COMMON if k in present]
    common = wide['partial_spearman'][listed].notna().all(axis=1) & wide['spearman'][listed].notna().all(axis=1)
    common = sorted(common[common].index)
    rows, per = [], []
    for k in present:
        for metric in ['spearman', 'partial_spearman']:
            col = wide[metric][k]
            on_common = col.reindex(common)
            defined = col.dropna()
            row = dict(slide=slide, rule=rule, level='receptor', condition=label(*k), metric=metric,
                       receptor_set='common', defined_in_condition=len(defined),
                       **boot(on_common.dropna().to_numpy()))
            row['missing_on_common'] = int(on_common.isna().sum())
            rows.append(row)
            per.append(on_common.rename(f'{label(*k)}|{metric}'))
    ro = r[r.receptor.isin(common)].groupby('receptor').receptor_only_spearman.agg(['min', 'max'])
    check(np.allclose(ro['min'], ro['max'], equal_nan=True), 'receptor-only Spearman differs between conditions')
    rows.append(dict(slide=slide, rule=rule, level='receptor', condition='receptor_only', metric='spearman',
                     receptor_set='common', defined_in_condition=int(ro['min'].notna().sum()),
                     missing_on_common=int(ro['min'].isna().sum()), **boot(ro['min'].dropna().to_numpy())))
    per.append(ro['min'].rename('receptor_only|spearman'))
    return pd.DataFrame(rows), pd.concat(per, axis=1).assign(slide=slide, rule=rule)


def pair_tables(m, lr, key, slide, rule):
    p = m[(m.level == 'pair') & (m.top_tf == key)].copy()
    kinds = {f'{a}->{b}': t for a, b, t in zip(lr.ligand, lr.receptor, lr.type)}
    p['lr_type'] = p.unit.map(kinds)
    check(p.lr_type.notna().all(), 'pair unit not in lr.csv')
    rows = []
    for k in ORDER:
        g = p[(p.condition == k[0]) & (p.r0_key == k[1])]
        if g.empty:
            continue
        for subset, gg in [('all', g), ('secreted', g[g.lr_type == 'diffusion']), ('contact', g[g.lr_type == 'contact'])]:
            for nan_rule in ['dropped', 'zero']:
                for metric in ['spearman', 'partial_spearman', 'receptor_only_spearman']:
                    v = gg[metric]
                    v = v.dropna() if nan_rule == 'dropped' else v.fillna(0.0)
                    rows.append(dict(slide=slide, rule=rule, level='pair', condition=label(*k), lr_subset=subset,
                                     nan=nan_rule, metric=metric, units_listed=len(gg), **boot(v.to_numpy())))
    return pd.DataFrame(rows)


def run(args):
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    launch = launch_info()
    rec, pair, per = [], [], []
    for slide, ev in [('posterior', Path(args.posterior)), ('anterior', Path(args.anterior))]:
        meta, m, lr = load(ev)
        for rule, key in [('main_top10', str(meta['primary_top_tf'])), ('commot_rule', 'commot')]:
            if not (m.top_tf == key).any():
                continue
            a, b = receptor_tables(m, key, slide, rule)
            rec.append(a)
            per.append(b)
            pair.append(pair_tables(m, lr, key, slide, rule))
    out.mkdir(parents=True)
    pd.concat(rec).to_csv(out/'receptor_level.csv', index=False)
    pd.concat(pair).to_csv(out/'pair_level.csv', index=False)
    pd.concat(per).to_csv(out/'receptor_values.csv')
    stage2.write_json(out/'manifest.json', dict(status='complete', launch=launch, posterior=args.posterior,
                      anterior=args.anterior, replicates=REPLICATES, seed=SEED, conditions=[label(*k) for k in ORDER],
                      common_set_conditions=[label(*k) for k in COMMON], note=__doc__))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--posterior', default='target_eval_stage2_post')
    p.add_argument('--anterior', default='target_eval_stage2_ant')
    p.add_argument('--output', default='reports/commot_eval_values')
    return p


if __name__ == '__main__':
    run(parser().parse_args())
