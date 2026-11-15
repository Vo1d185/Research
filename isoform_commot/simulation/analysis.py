"""Section 9 of reports/pde_r0_plan.md, from the per-pair metrics of `simulation.benchmark evaluate`.

Settings: main (Spearman, steady-state truth, eps/rho 0.1/10, generated scenarios) and the secondary
analyses t50 (Spearman, t = 50 truth), eps0.2_rho1 (Spearman, steady), rmse (RMSE, steady, 0.1/10).
Aggregation: mean over pairs (scenario) -> mean over repeats (structure) -> equal-weight mean over structures.
r0*(lambda) = best of the 7 r0/lambda values (Spearman: max; RMSE: min); exact ties -> smaller r0.
H holds iff r0* (absolute) is non-decreasing in lambda and r0*(largest) > r0*(smallest).
Runs failing the hard convergence requirements (stage-2 8.7) are excluded and counted.
Bootstrap (descriptive): repeats resampled with replacement within each structure, 1000 times, default_rng(1).
"""
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import kendalltau

from simulation.benchmark import DIS_FIXED, DIS_OVER_LAMBDA, H, KS, R0_OVER_LAMBDA, fmt, lam

GRID = [f'r0_{fmt(g)}' for g in R0_OVER_LAMBDA]
SETTINGS = {
    'main': dict(metric='spearman', truth='steady', eps=0.1, rho=10.0),
    't50': dict(metric='spearman', truth='t50', eps=0.1, rho=10.0),
    'eps0.2_rho1': dict(metric='spearman', truth='steady', eps=0.2, rho=1.0),
    'rmse': dict(metric='rmse', truth='steady', eps=0.1, rho=10.0),
}
REPLICATES, SEED = 1000, 1


def converged(m):
    passed = m.check_passed.astype(str).str.lower().eq('true')
    return passed & (pd.to_numeric(m.kkt_tol, errors='coerce') <= 1e-7) & (pd.to_numeric(m.forbidden_mass) == 0)


def label(m):
    """Analysis label per row: inv, r0_<g>, cd_<m>lambda (m * lambda_k), cd_150, ref; other cd rows dropped."""
    parts = []
    for k in KS:
        sub = m[m.k == k]
        o = sub[sub.family.isin(['inv', 'r0'])].copy()
        o['cond'] = o.condition
        parts.append(o)
        cd = sub[sub.family == 'cd']
        for mult in DIS_OVER_LAMBDA:
            if mult * lam(k) >= H:
                x = cd[np.isclose(cd.dis_thr, mult * lam(k), rtol=0, atol=1e-6)].copy()
                x['cond'] = f'cd_{fmt(mult)}lambda'
                parts.append(x)
        x = cd[np.isclose(cd.dis_thr, DIS_FIXED)].copy()
        x['cond'] = f'cd_{fmt(DIS_FIXED)}'
        parts.append(x)
        x = sub[sub.family == 'ref'].copy()
        x['cond'] = 'ref'
        parts.append(x)
    return pd.concat(parts, ignore_index=True)


def tensor(m, metric, conds):
    """Scenario means (mean over pairs) as an array [structure, repeat, k, cond] (NaN where missing)."""
    sm = m.groupby(['structure', 'repeat', 'k', 'cond'])[metric].mean()
    structures, repeats = sorted(m.structure.unique()), sorted(m.repeat.unique())
    idx = pd.MultiIndex.from_product([structures, repeats, list(KS), conds])
    return sm.reindex(idx).to_numpy().reshape(len(structures), len(repeats), len(KS), len(conds)), structures


def level(t, pick=None):
    """lambda-level means [k, cond] from a scenario tensor; pick: repeat indices per structure."""
    if pick is not None:
        t = np.take_along_axis(t, pick[:, :, None, None], axis=1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)          # all-NaN slices stay NaN
        return np.nanmean(np.nanmean(t, axis=1), axis=0)


def best(values, metric):
    """Index of the best grid value (ascending grid): first max (Spearman) or first min (RMSE)."""
    return int(np.nanargmax(values) if metric == 'spearman' else np.nanargmin(values))


def exponents(choice):
    """log2 of the absolute r0* in units of lambda_A: log2(g*) + k (grid and lambda both powers of 2)."""
    return [int(round(np.log2(R0_OVER_LAMBDA[c]))) + k for c, k in zip(choice, KS)]


def holds(choice):
    n = exponents(choice)
    return bool(all(b >= a for a, b in zip(n, n[1:])) and n[-1] > n[0])


def analyse_setting(m, name, setting, rng_seed=SEED):
    metric = setting['metric']
    sel = m[(m.source == 'generated') & (m.truth == setting['truth']) & np.isclose(m.eps, setting['eps']) &
            np.isclose(m.rho, setting['rho'])]
    t, structures = tensor(sel, metric, GRID)
    lv = level(t)
    choice = [best(lv[i], metric) for i in range(len(KS))]
    rows = []
    for i, k in enumerate(KS):
        g = R0_OVER_LAMBDA[choice[i]]
        rows.append(dict(setting=name, k=k, lambda_over_h=lam(k) / H, r0_over_lambda=g, r0_over_h=g * lam(k) / H,
                         at_boundary=g in (R0_OVER_LAMBDA[0], R0_OVER_LAMBDA[-1]),
                         **{c: lv[i, j] for j, c in enumerate(GRID)}))
    selection = pd.DataFrame(rows)
    verdict = holds(choice)
    rng = np.random.default_rng(rng_seed)
    reps = t.shape[1]
    boot_choice, boot_level = [], []
    for _ in range(REPLICATES):
        pick = rng.integers(0, reps, size=(t.shape[0], reps))
        bl = level(t, pick)
        boot_level.append(bl)
        boot_choice.append([best(bl[i], metric) for i in range(len(KS))])
    boot_choice, boot_level = np.array(boot_choice), np.array(boot_level)
    boot = []
    for i, k in enumerate(KS):
        g = np.array(R0_OVER_LAMBDA)[boot_choice[:, i]]
        boot.append(dict(setting=name, k=k, ci_low=float(np.quantile(g, 0.025, method='inverted_cdf')),
                         ci_high=float(np.quantile(g, 0.975, method='inverted_cdf')),
                         **{f'share_{fmt(x)}': float(np.mean(g == x)) for x in R0_OVER_LAMBDA}))
    boot = pd.DataFrame(boot)
    se = pd.DataFrame([dict(setting=name, k=k, cond=c, mean=lv[i, j], se=float(np.nanstd(boot_level[:, i, j], ddof=1)))
                       for i, k in enumerate(KS) for j, c in enumerate(GRID)])
    per_structure = []
    xs, ys = [], []
    for s_i, s in enumerate(structures):
        st = np.nanmean(t[s_i], axis=0)
        ch = [best(st[i], metric) for i in range(len(KS))]
        n = exponents(ch)
        xs += list(KS); ys += n
        per_structure.append(dict(setting=name, structure=s, holds=holds(ch),
                                  **{f'r0_over_lambda_k{k:+d}': R0_OVER_LAMBDA[c] for k, c in zip(KS, ch)}))
    tau = kendalltau(xs, ys)
    summary = dict(setting=name, metric=metric, truth=setting['truth'], eps=setting['eps'], rho=setting['rho'],
                   hypothesis_holds=verdict, r0_over_lambda={str(k): R0_OVER_LAMBDA[c] for k, c in zip(KS, choice)},
                   bootstrap_share_holding=float(np.mean([holds(c) for c in boot_choice])),
                   structures_holding=int(sum(r['holds'] for r in per_structure)), structures=len(per_structure),
                   kendall_tau=float(tau.statistic), kendall_p=float(tau.pvalue),
                   missing_scenario_values=int(np.isnan(t).sum()),
                   spearman_nan_pairs=int(sel.spearman_nan.astype(str).str.lower().eq('true').sum()))
    return summary, selection, boot, se, pd.DataFrame(per_structure)


def comparisons(m, setting=SETTINGS['main'], rng_seed=SEED):
    """Descriptive: COMMOT dis_thr*, cost-only differences against Cd(4 lambda), method table (plan section 9)."""
    metric = setting['metric']
    sel = m[(m.source == 'generated') & (m.truth == setting['truth']) & np.isclose(m.eps, setting['eps']) &
            np.isclose(m.rho, setting['rho'])]
    cd = [f'cd_{fmt(x)}lambda' for x in DIS_OVER_LAMBDA]
    conds = ['inv'] + GRID + cd + [f'cd_{fmt(DIS_FIXED)}']
    t, _ = tensor(sel, metric, conds)
    lv = level(t)
    ci = {c: j for j, c in enumerate(conds)}
    rng = np.random.default_rng(rng_seed)
    boots = []
    for _ in range(REPLICATES):
        pick = rng.integers(0, t.shape[1], size=(t.shape[0], t.shape[1]))
        boots.append(np.take_along_axis(t, pick[:, :, None, None], axis=1))
    rows, diffs = [], []
    for i, k in enumerate(KS):
        g = best(lv[i, [ci[c] for c in GRID]], metric)
        avail = [c for c in cd if not np.isnan(lv[i, ci[c]])]
        dstar = avail[best(np.array([lv[i, ci[c]] for c in avail]), metric)] if avail else None
        table = {'O-r0 at r0* (selected)': GRID[g], 'O-r0 at r0/lambda = 1': 'r0_1', 'O-inv': 'inv',
                 'Cd at dis_thr* (selected)': dstar, 'Cd at 4 lambda': 'cd_4lambda', 'Cd at 150': f'cd_{fmt(DIS_FIXED)}'}
        for label_, c in table.items():
            if c is None:
                continue
            j = ci[c]
            se = np.nanstd([np.nanmean(np.nanmean(b[:, :, i, j], axis=1)) for b in boots], ddof=1)
            rows.append(dict(k=k, setting=label_, condition=c, mean=lv[i, j], bootstrap_se=se,
                             selected='selected' in label_))
        for name, a in [('O-r0 at r0* - Cd 4 lambda', GRID[g]), ('O-r0 at r0/lambda = 1 - Cd 4 lambda', 'r0_1'),
                        ('O-inv - Cd 4 lambda', 'inv')]:
            d = t[:, :, i, ci[a]] - t[:, :, i, ci['cd_4lambda']]
            est = np.nanmean(np.nanmean(d, axis=1))
            bd = [np.nanmean(np.nanmean(np.take_along_axis(d, rng.integers(0, d.shape[1], size=d.shape), axis=1), axis=1))
                  for _ in range(REPLICATES)]
            diffs.append(dict(k=k, comparison=name, metric=metric, mean_difference=est,
                              ci_low=float(np.quantile(bd, 0.025)), ci_high=float(np.quantile(bd, 0.975))))
        if dstar:
            rows.append(dict(k=k, setting='COMMOT dis_thr*/lambda', condition=dstar, mean=np.nan, bootstrap_se=np.nan,
                             selected=True))
    return pd.DataFrame(rows), pd.DataFrame(diffs)


def archive_comparison(m):
    """Plan 4.3: lambda_A, eps/rho 0.1/10: archived (one per structure) vs generated scenarios, both truths."""
    rows = []
    conds = ['inv'] + GRID + [f'cd_{fmt(x)}lambda' for x in DIS_OVER_LAMBDA] + [f'cd_{fmt(DIS_FIXED)}']
    for truth in ('steady', 't50'):
        for source in ('archive', 'generated'):
            sel = m[(m.source == source) & (m.truth == truth) & (m.k == 0) & np.isclose(m.eps, 0.1) & np.isclose(m.rho, 10)]
            for metric in ('spearman', 'rmse'):
                st = sel.groupby(['structure', 'repeat', 'cond'])[metric].mean().groupby(['structure', 'cond']).mean()
                lv = st.groupby('cond').mean().reindex(conds)
                g = GRID[best(lv[GRID].to_numpy(), metric)]
                rows.append(dict(truth=truth, source=source, metric=metric, r0_star=g, **lv.to_dict()))
    return pd.DataFrame(rows)


def run(metrics_path, output, runs_dir=None):
    m = pd.read_csv(metrics_path)
    m = label(m)
    ok = converged(m)
    excluded = m[~ok].drop_duplicates(['scenario', 'condition', 'eps', 'rho', 'nitermax'])
    m = m[ok]
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    excluded.to_csv(out/'excluded_runs.csv', index=False)
    summaries, parts = [], {k: [] for k in ['selection', 'bootstrap', 'curves', 'per_structure']}
    for name, setting in SETTINGS.items():
        s, sel, boot, se, per = analyse_setting(m, name, setting)
        summaries.append(s)
        for key, v in zip(['selection', 'bootstrap', 'curves', 'per_structure'], [sel, boot, se, per]):
            parts[key].append(v)
    for key, v in parts.items():
        pd.concat(v, ignore_index=True).to_csv(out/f'{key}.csv', index=False)
    table, diffs = comparisons(m)
    table.to_csv(out/'methods.csv', index=False)
    diffs.to_csv(out/'cost_only.csv', index=False)
    archive_comparison(m).to_csv(out/'archive_comparison.csv', index=False)
    (out/'analysis.json').write_text(json.dumps(dict(summaries=summaries, excluded_runs=len(excluded)), indent=2,
                                                default=str))
    return summaries
