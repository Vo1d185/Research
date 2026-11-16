"""Stage-2 plan section 9: is the r0 kernel useful? Tables and figures from one target evaluation.

  posterior (selection, exploration):
    python r0_validation.py --eval target_eval_stage2_post --output reports/stage2_r0val_post
  anterior (confirmatory test, r0* fixed by the posterior):
    python r0_validation.py --eval target_eval_stage2_ant --output reports/stage2_r0val_ant \
        --posterior reports/stage2_r0val_post/validation.json

Conditions (stage2.condition_cases): A (r0 kernel, four r0), F (flat kernel), I (inverse kernel),
B (COMMOT cost, pure collective), P1-P3 (posterior only: positions permuted, r0 = 200 um).
Metric: receptor-level partial Spearman (receptor expression controlled), top-10 TFs, mean target expression.
Intervals: bootstrap over receptors, 1000 resamples, numpy.random.default_rng(1), 2.5% and 97.5% quantiles.
A group of comparisons uses the receptors defined in every condition of the group (9.3).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

import stage2
from research import check, launch_info

REPLICATES, SEED = 1000, 1
REFERENCES = ('F', 'I', 'B')
SHUFFLED = ('P1', 'P2', 'P3')


def clean(obj):
    """JSON-safe copy: numpy scalars -> Python, NaN -> None."""
    if isinstance(obj, dict):
        return {str(k): clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return None if np.isnan(obj) else float(obj)
    return obj


def bootstrap_mean(values, seed=SEED):
    """Mean and percentile interval of the mean, resampling receptors with replacement."""
    x = np.asarray(values, float)
    if len(x) == 0:
        return dict(mean=np.nan, ci_low=np.nan, ci_high=np.nan, n=0)
    idx = np.random.default_rng(seed).integers(0, len(x), size=(REPLICATES, len(x)))
    means = x[idx].mean(axis=1)
    return dict(mean=float(x.mean()), ci_low=float(np.quantile(means, 0.025)), ci_high=float(np.quantile(means, 0.975)),
                n=len(x))


def bootstrap_group_difference(x, y, seed=SEED):
    """mean(x) - mean(y) for two receptor groups, each resampled on its own (plan 9.5 d)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    if len(x) == 0 or len(y) == 0:
        return dict(difference=np.nan, ci_low=np.nan, ci_high=np.nan)
    rng = np.random.default_rng(seed)
    bx = x[rng.integers(0, len(x), size=(REPLICATES, len(x)))].mean(axis=1)
    by = y[rng.integers(0, len(y), size=(REPLICATES, len(y)))].mean(axis=1)
    d = bx - by
    return dict(difference=float(x.mean() - y.mean()), ci_low=float(np.quantile(d, 0.025)),
                ci_high=float(np.quantile(d, 0.975)))


def columns(wide, r0s):
    """Series per analysis label: A at each r0, F, I, B, P1-P3 (those present)."""
    out = {}
    for r0 in r0s:
        if ('A', r0) in wide.columns:
            out[f'A_{r0:g}'] = wide[('A', r0)]
    for code in REFERENCES + SHUFFLED:
        if code in wide.columns.get_level_values(0):
            sub = wide[code]
            check(sub.shape[1] == 1, f'{code}: expected one case')
            out[code] = sub.iloc[:, 0]
    return pd.DataFrame(out)


def curve_table(table, r0s, receptor_only):
    """Mean partial Spearman with interval on the receptors defined in A (all r0), F, I and B."""
    need = [f'A_{r0:g}' for r0 in r0s] + [c for c in REFERENCES if c in table]
    common = table[need].dropna().index
    rows = []
    for name in need:
        kind, r0 = ('A', float(name[2:])) if name.startswith('A_') else (name, np.nan)
        rows.append(dict(condition=kind, r0_um=r0, **bootstrap_mean(table.loc[common, name])))
    rows.append(dict(condition='receptor_only', r0_um=np.nan, **bootstrap_mean(receptor_only.reindex(common))))
    return pd.DataFrame(rows), common


def difference_tables(table, r0_star, classes):
    """Paired differences A(r0*) - F / I / B on the receptors defined in all four; by receptor class."""
    from scipy.stats import binomtest
    a = f'A_{r0_star:g}'
    refs = [c for c in REFERENCES if c in table]
    common = table[[a] + refs].dropna().index
    cls = pd.Series(common.map(classes), index=common)
    summary, receptors = [], []
    for ref in refs:
        d = (table.loc[common, a] - table.loc[common, ref]).rename('difference')
        receptors.append(pd.DataFrame(dict(reference=ref, receptor=common, receptor_class=cls.to_numpy(),
                                           a=table.loc[common, a].to_numpy(), other=table.loc[common, ref].to_numpy(),
                                           difference=d.to_numpy())))
        for group in ['all', 'Secreted', 'Contact', 'both']:
            g = d if group == 'all' else d[cls == group]
            pos, neg = int((g > 0).sum()), int((g < 0).sum())
            summary.append(dict(comparison=f'A(r0*) - {ref}', reference=ref, receptor_class=group,
                                **bootstrap_mean(g), a_higher=pos, reference_higher=neg, ties=int((g == 0).sum()),
                                sign_test_p=float(binomtest(pos, pos + neg, 0.5).pvalue) if pos + neg else np.nan))
    return pd.DataFrame(summary), pd.concat(receptors, ignore_index=True), common


def criteria(summary, receptors):
    """Plan 9.5 on one slide (the verdict is used for the anterior slide only)."""
    row = lambda ref, group='all': summary[(summary.reference == ref) & (summary.receptor_class == group)].iloc[0]
    excludes_zero = lambda r: bool(r['mean'] > 0 and r['ci_low'] > 0)
    out = {}
    f = row('F')
    out['a'] = dict(met=excludes_zero(f), mean_difference=f['mean'], ci_low=f['ci_low'], ci_high=f['ci_high'], n=int(f['n']))
    parts = {ref: row(ref) for ref in ('I', 'B')}
    out['b'] = dict(met=all(excludes_zero(r) for r in parts.values()),
                    **{f'{ref}_{k}': (bool(excludes_zero(r)) if k == 'met' else float(r[k]))
                       for ref, r in parts.items() for k in ('met', 'mean', 'ci_low', 'ci_high')})
    out['c'] = dict(met=bool(f['a_higher'] > f['reference_higher']), a_higher=int(f['a_higher']),
                    f_higher=int(f['reference_higher']), ties=int(f['ties']), sign_test_p=f['sign_test_p'])
    d = receptors[receptors.reference == 'F']
    sec, con = d[d.receptor_class == 'Secreted'].difference, d[d.receptor_class == 'Contact'].difference
    gd = bootstrap_group_difference(sec, con)
    out['d'] = dict(met=bool(len(sec) and len(con) and sec.mean() > con.mean()), secreted_mean=float(sec.mean()),
                    contact_mean=float(con.mean()), secreted_n=len(sec), contact_n=len(con),
                    both_mean=float(d[d.receptor_class == 'both'].difference.mean()),
                    both_n=int((d.receptor_class == 'both').sum()), **gd)
    met = [out[k]['met'] for k in 'abcd']
    out['verdict'] = 'useful' if all(met) else 'no effect visible under this evaluation' if not any(met) else \
        'partly met: ' + ', '.join(k for k, m in zip('abcd', met) if m)
    return out


def negative_control(table):
    """Plan 9.6: A(r0 = 200 um) minus the mean over the three permuted runs, per receptor."""
    a = f'A_{stage2.SHUFFLE_R0_UM:g}'
    present = [c for c in SHUFFLED if c in table]
    if a not in table or len(present) != len(SHUFFLED):
        return None, None
    common = table[[a] + present].dropna().index
    rows = [dict(comparison='A(200) - mean of P1-P3', **bootstrap_mean(table.loc[common, a] - table.loc[common, present].mean(axis=1)))]
    for c in present:
        rows.append(dict(comparison=f'A(200) - {c}', **bootstrap_mean(table.loc[common, a] - table.loc[common, c])))
    rows.append(dict(comparison='A(200) mean partial Spearman', **bootstrap_mean(table.loc[common, a])))
    for c in present:
        rows.append(dict(comparison=f'{c} mean partial Spearman', **bootstrap_mean(table.loc[common, c])))
    out = pd.DataFrame(rows)
    main = out.iloc[0]
    insensitive = bool(main.ci_low <= 0 <= main.ci_high)
    return out, dict(insensitive=insensitive, mean_difference=float(main['mean']), ci_low=float(main.ci_low),
                     ci_high=float(main.ci_high), receptors=int(main.n))


def figure_curve(curve, r0_star, title, path):
    """Left: A over r0 with intervals, reference means as thin lines. Right: the reference conditions with intervals."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from analyze import PALETTE, style_axis
    colors = dict(A=PALETTE['series1'], B=PALETTE['series2'], F=PALETTE['series3'], I='#eda100',
                  receptor_only=PALETTE['muted'])
    labels = dict(A='A: r0 kernel', B='B: COMMOT cost', F='F: flat kernel', I='I: inverse kernel',
                  receptor_only='receptor expression only (plain Spearman)')
    short = dict(B='B\nCOMMOT', F='F\nflat', I='I\ninverse', receptor_only='receptor\nonly')
    styles = dict(B='--', F='-.', I=':', receptor_only=(0, (1, 3)))
    marks = dict(B='s', F='D', I='^', receptor_only='v')
    fig, (ax, rx) = plt.subplots(1, 2, figsize=(8.4, 4.2), constrained_layout=True, sharey=True,
                                 gridspec_kw=dict(width_ratios=[3, 2]))
    a = curve[curve.condition == 'A'].sort_values('r0_um')
    x = np.log2(a.r0_um.to_numpy())
    span = (x.min() - 0.35, x.max() + 0.35)
    refs = [c for c in ['F', 'I', 'B', 'receptor_only'] if (curve.condition == c).any()]
    handles = []
    for k, code in enumerate(refs):
        row = curve[curve.condition == code].iloc[0]
        ax.plot(span, [row['mean']]*2, color=colors[code], ls=styles[code], lw=1.5)
        h = rx.errorbar([k], [row['mean']], yerr=[[row['mean'] - row.ci_low], [row.ci_high - row['mean']]], color=colors[code],
                        marker=marks[code], ms=8, mec=PALETTE['surface'], mew=1.5, lw=2, capsize=3, label=labels[code])
        handles.append(h)
    h = ax.errorbar(x, a['mean'], yerr=[a['mean'] - a.ci_low, a.ci_high - a['mean']], color=colors['A'], lw=2, marker='o',
                    ms=8, mec=PALETTE['surface'], mew=1.5, capsize=3, label=labels['A'], zorder=3)
    ax.set_xlim(*span)
    ax.set_xticks(x, [f'{v:g}' for v in a.r0_um])
    ax.set_xlabel('r0 (um)', fontsize=9, color=PALETTE['text_secondary'])
    ax.set_ylabel(f'Mean partial Spearman ({int(a.n.iloc[0])} receptors)', fontsize=9, color=PALETTE['text_secondary'])
    ax.set_title(f'A: r0 kernel, r0* = {r0_star:g} um (lines: reference means)', fontsize=9, loc='left')
    rx.set_xlim(-0.6, len(refs) - 0.4)
    rx.set_xticks(range(len(refs)), [short[c] for c in refs])
    rx.set_title('Reference conditions', fontsize=9, loc='left')
    for axis in (ax, rx):
        style_axis(axis)
    fig.suptitle(title + '\n(bars: 95% bootstrap intervals over receptors)', fontsize=10, x=0.01, ha='left')
    fig.legend(handles=[h] + handles, fontsize=8, frameon=False, loc='outside lower center', ncol=3)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def figure_differences(receptors, summary, title, path):
    """Per-receptor differences A(r0*) - F and A(r0*) - B, Secreted and Contact receptors in separate panels."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from analyze import PALETTE, style_axis
    refs, groups = ['F', 'B'], ['Secreted', 'Contact']
    fig, axes = plt.subplots(len(refs), len(groups), figsize=(8.4, 5.6), constrained_layout=True, sharey='row',
                             squeeze=False)
    for i, ref in enumerate(refs):
        for j, group in enumerate(groups):
            ax = axes[i, j]
            d = receptors[(receptors.reference == ref) & (receptors.receptor_class == group)].difference.sort_values()
            s = summary[(summary.reference == ref) & (summary.receptor_class == group)]
            ax.axhline(0, color=PALETTE['axis'], lw=1)
            if len(d):
                r = s.iloc[0]
                ax.axhspan(r.ci_low, r.ci_high, color=PALETTE['series1'], alpha=.12, lw=0)
                ax.axhline(r['mean'], color=PALETTE['series1'], lw=2)
                ax.scatter(np.arange(len(d)), d, s=22, color=PALETTE['series1'], edgecolor=PALETTE['surface'], lw=.8, zorder=3)
                ax.text(.02, .96, f"mean {r['mean']:+.4f}  (95% CI {r.ci_low:+.4f} to {r.ci_high:+.4f})\n"
                        f"{int(r.a_higher)} higher, {int(r.reference_higher)} lower, n = {int(r.n)}",
                        transform=ax.transAxes, va='top', fontsize=8, color=PALETTE['text_secondary'])
            ax.set_title(f'A(r0*) - {ref}: {group} receptors', fontsize=9, loc='left')
            ax.set_xlabel('receptors, sorted by difference', fontsize=8, color=PALETTE['text_secondary'])
            ax.set_xticks([])
            if j == 0:
                ax.set_ylabel('Difference in partial Spearman', fontsize=8, color=PALETTE['text_secondary'])
            style_axis(ax)
    fig.suptitle(title, fontsize=10, x=0.01, ha='left')
    fig.savefig(path, dpi=180)
    plt.close(fig)


def analyse(ev, meta, cases, classes, rule, r0_star=None):
    m, wide, evaluable = stage2.metric_table(ev, meta, cases, rule)
    selection, own_r0 = stage2.select_r0(wide)
    r0s = sorted(float(v) for v in wide['A'].columns)
    table = columns(wide, r0s)
    receptor_only = m[m.condition == 'A'].groupby('receptor').receptor_only_spearman.first()
    star = own_r0 if r0_star is None else r0_star
    curve, curve_set = curve_table(table, r0s, receptor_only)
    summary, receptors, diff_set = difference_tables(table, star, classes)
    control_table, control = negative_control(table)
    return dict(selection=selection, own_r0=own_r0, r0_star=star, curve=curve, summary=summary, receptors=receptors,
                control_table=control_table, control=control, criteria=criteria(summary, receptors),
                receptors_evaluable=len(evaluable), receptors_curve=len(curve_set), receptors_differences=len(diff_set))


def run(args):
    launch = launch_info()
    ev, out = Path(args.eval), Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    meta = json.loads((ev/'eval_manifest.json').read_text())
    check(meta.get('status') == 'complete', 'Incomplete target evaluation')
    cases = pd.concat([stage2.condition_cases(r) for r in meta['results']], ignore_index=True)
    cases = cases[cases.condition.isin(('A', 'B') + stage2.VALIDATION)].reset_index(drop=True)
    check(not cases.duplicated(['condition', 'r0_um'], keep=False).any(), 'A condition/r0 appears twice')
    check({'A', 'B', 'F', 'I'} <= set(cases.condition), 'Need conditions A, B, F and I')
    for case in cases.to_dict('records'):                       # plan 8.7 / 9.2 hard requirements
        case['forbidden_mass'] = stage2.forbidden_mass(case)
        check(stage2.converged(case) is True, f"{case['results']}/{case['case']}: fails the convergence requirements")
    classes = stage2.receptor_classes(pd.read_csv(Path(meta['input'])/'lr.csv', dtype=str))
    posterior = None
    if args.posterior:
        posterior = json.loads(Path(args.posterior).read_text())
        check(posterior.get('status') == 'complete' and posterior.get('role') == 'selection', 'Need the posterior validation')
    role = 'confirmation' if posterior else 'selection'
    r0_fixed = float(posterior['r0_star_um']) if posterior else None
    main = analyse(ev, meta, cases, classes, 'ours', r0_fixed)
    r0_star = main['r0_star']
    out.mkdir(parents=True)
    slide = 'anterior' if posterior else 'posterior'
    main['selection'].to_csv(out/'r0_selection.csv', index=False)
    main['curve'].to_csv(out/'curve.csv', index=False)
    main['summary'].to_csv(out/'differences.csv', index=False)
    main['receptors'].to_csv(out/'differences_receptors.csv', index=False)
    if main['control_table'] is not None:
        main['control_table'].to_csv(out/'negative_control.csv', index=False)
    figure_curve(main['curve'], r0_star, f'{slide.capitalize()} slide: mean partial Spearman over r0, with reference conditions',
                 out/'curve.png')
    figure_differences(main['receptors'], main['summary'],
                       f'{slide.capitalize()} slide: per-receptor differences at r0* = {r0_star:g} um', out/'differences.png')
    sets = pd.read_csv(ev/'target_sets.csv', dtype={'top_tf': str})
    commot = None
    if 'rule' in sets and (sets.rule == 'commot').any():        # plan 8.8 / 9.3: sensitivity, same r0*
        c = analyse(ev, meta, cases, classes, 'commot', r0_star)
        c['curve'].to_csv(out/'curve_commot_rule.csv', index=False)
        c['summary'].to_csv(out/'differences_commot_rule.csv', index=False)
        if c['control_table'] is not None:
            c['control_table'].to_csv(out/'negative_control_commot_rule.csv', index=False)
        commot = dict(criteria=c['criteria'], negative_control=c['control'], receptors_curve=c['receptors_curve'],
                      receptors_differences=c['receptors_differences'], rule_selected_r0_um=c['own_r0'])
    stage2.write_json(out/'validation.json', clean(dict(status='complete', launch=launch, role=role, slide=slide, eval=str(ev),
        input=meta['input'], results=sorted(cases.results.unique()), top_tf=meta['primary_top_tf'],
        target_score=meta['target_score'][0], r0_star_um=r0_star, rule_selected_r0_um=main['own_r0'],
        r0_rule='plan section 8.5 (posterior); fixed on the anterior', replicates=REPLICATES, seed=SEED,
        receptors_evaluable=main['receptors_evaluable'], receptors_curve=main['receptors_curve'],
        receptors_differences=main['receptors_differences'], receptor_classes=classes.value_counts().to_dict(),
        criteria=main['criteria'], criteria_role='confirmatory (plan 9.5)' if posterior else 'exploratory: r0* selected on this slide',
        negative_control=main['control'], commot_rule=commot,
        eval_sha256={f: stage2.digest(ev/f) for f in ['case_metrics.csv', 'target_sets.csv']})))
    print(main['curve'].to_string(index=False))
    print(main['summary'].to_string(index=False))
    print(json.dumps(clean(dict(r0_star_um=r0_star, criteria=main['criteria'], negative_control=main['control'])), indent=2))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--eval', required=True, help='evaluate_targets.py output with A, B, F, I (posterior: P1-P3)')
    p.add_argument('--output', required=True)
    p.add_argument('--posterior', help='anterior only: validation.json of the posterior (fixes r0*)')
    return p


if __name__ == '__main__':
    run(parser().parse_args())
