"""Analyze saved transport plans; never rerun OT or change the cost."""
import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from research import check, numeric_table


# Validated reference palette (first three categorical slots pass all-pairs CVD checks);
# slot 3 is below 3:1 contrast, so it is always paired with a distinct marker + legend.
PALETTE = dict(series1='#2a78d6', series2='#eb6834', series3='#1baf7a', muted='#898781',
               text='#0b0b0b', text_secondary='#52514e', grid='#e1e0d9', axis='#c3c2b7',
               surface='#fcfcfb')
TYPE_COLORS = {'diffusion': PALETTE['series1'], 'contact': PALETTE['series2']}


def style_axis(ax):
    ax.set_facecolor(PALETTE['surface'])
    ax.grid(color=PALETTE['grid'], lw=.6)
    ax.set_axisbelow(True)
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    for side in ['left', 'bottom']:
        ax.spines[side].set_color(PALETTE['axis'])
    ax.tick_params(colors=PALETTE['text_secondary'], labelsize=8)
    ax.title.set_color(PALETTE['text'])


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def decompose(base, scored, blocks):
    """Split allocation change into between-LR-pair shares and within-pair spatial patterns.

    blocks: (row_slice, col_slice) per LR pair. between_TV compares each pair's share of
    total transport; within_TV is the share-weighted TV of each pair's own normalized
    spatial plan. A pure rescaling of one pair changes between_TV only.
    """
    x, y = float(base.sum()), float(scored.sum())
    if not x or not y:
        return dict(between_pair_TV=None, within_pair_TV=None)
    mb = np.array([base[r, c].sum() for r, c in blocks])
    ms = np.array([scored[r, c].sum() for r, c in blocks])
    qb, qs = mb/x, ms/y
    within = 0.
    for (r, c), b_mass, s_mass, wb, ws in zip(blocks, mb, ms, qb, qs):
        if b_mass > 0 and s_mass > 0:
            within += 0.5*(wb + ws) * 0.5*np.abs(scored[r, c]/s_mass - base[r, c]/b_mass).sum()
    return dict(between_pair_TV=float(0.5*np.abs(qs - qb).sum()), within_pair_TV=float(within))


def paired_metrics(base, scored):
    """Both arrays must use exactly the same entities/spots and distance design."""
    x, y = float(base.sum()), float(scored.sum())
    return dict(
        baseline_total_T=x, ppi_total_T=y,
        total_change_percent=100 * (y/x - 1) if x else None,
        paired_relative_L1=float(np.abs(scored-base).sum()/x) if x else None,
        allocation_TV=float(np.abs(scored/y-base/x).sum()/2) if x and y else None,
    )


def type_and_pair_metrics(T, C, spots, ligands, receptors, pairs, D):
    n = len(spots)
    check(T.shape == C.shape == (len(ligands)*n, len(receptors)*n), 'Plan shape mismatch')
    check(np.isfinite(T).all() and (T >= 0).all(), 'Invalid transport')
    check((T[~np.isfinite(C)] <= 1e-12).all(), 'Transport on forbidden connections')
    li, ri = {v:i for i,v in enumerate(ligands)}, {v:i for i,v in enumerate(receptors)}
    totals = {kind: [0., 0., 0] for kind in ['diffusion', 'contact']}
    lr_rows = []
    for row in pairs.itertuples(index=False):
        check(row.type in totals and row.ligand in li and row.receptor in ri,
              'Pair type or entity does not match saved plan')
        i, j = li[row.ligand], ri[row.receptor]
        block = T[i*n:(i+1)*n, j*n:(j+1)*n]
        mass = float(block.sum())
        distance_sum = float((block*D).sum())
        totals[row.type][0] += mass
        totals[row.type][1] += distance_sum
        totals[row.type][2] += int(np.isfinite(C[i*n:(i+1)*n,j*n:(j+1)*n]).sum())
        lr_rows.append(dict(ligand=row.ligand, receptor=row.receptor,
                            type=row.type, total_T=mass,
                            mean_distance=distance_sum/mass if mass else None))
    check(np.isclose(sum(v[0] for v in totals.values()), T.sum(), atol=1e-8),
          'LR/type partition does not account for all transport')
    rows = [dict(type=k, total_T=v[0], mean_distance=v[1]/v[0] if v[0] else None,
                 fraction_total_T=v[0]/T.sum() if T.sum() else None,
                 eligible_edges=v[2]) for k,v in totals.items()]
    return rows, lr_rows


def plot_results(types, paired, changes, out, unit, level):
    colors = TYPE_COLORS
    families = [f for f in ['r0_raw', 'r0_shape'] if (types.family == f).any()]
    fig, axes = plt.subplots(max(len(families), 1), 2, figsize=(11, 3.5*max(len(families), 1)),
                             constrained_layout=True, squeeze=False)
    if not families:
        for ax in axes.ravel():
            ax.text(.5, .5, 'No r0 design in these results', transform=ax.transAxes, ha='center')
            ax.set_axis_off()
    for row, family in enumerate(families):
        for col, metric in enumerate(['total_T', 'mean_distance']):
            ax = axes[row,col]
            for kind, color in colors.items():
                data = types[(types.family == family) & (types.score_mode == 'no_ppi')
                             & (types.type == kind)].sort_values('r0_h')
                ax.plot(data.r0_h, data[metric], 'o-', label=kind, color=color)
                base = types[(types.family == 'inverse') & (types.score_mode == 'no_ppi')
                             & (types.type == kind)]
                if len(base) and pd.notna(base.iloc[0][metric]):
                    ax.axhline(base.iloc[0][metric], color=color, ls='--', alpha=.7,
                               label=f'{kind}: inverse baseline')
            ax.set_title(f'{family} / no PPI')
            ax.set_xlabel('r0 / median nearest-neighbor spacing')
            ax.set_ylabel('Transport amount' if col == 0 else f'Mean distance ({unit})')
            ax.legend(fontsize=8, frameon=False)
            style_axis(ax)
            ax.margins(y=.12)
    fig.suptitle(f'Distance sensitivity by interaction type | {level}')
    fig.savefig(out/'distance_by_type.png', dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(11, 7), constrained_layout=True)
    if len(paired):
        labels = [f'{r.family}\n{r.r0_h:g} h' if pd.notna(r.r0_h) else r.family
                  for r in paired.itertuples()]
        for ax, metric, label in zip(axes, ['total_change_percent', 'allocation_TV'],
                                     ['Change in total T (%)', 'Allocation TV (0 to 1)']):
            ax.bar(np.arange(len(paired)), paired[metric], color=PALETTE['series1'], width=.6)
            ax.set_xticks(np.arange(len(paired)), labels, fontsize=8)
            ax.set_ylabel(label)
            ax.axhline(0, color=PALETTE['axis'], lw=.8)
            style_axis(ax)
            if metric == 'allocation_TV':
                ax.set_ylim(0, 1)
    else:
        for ax in axes:
            ax.text(.5, .5, 'No PPI input: matched PPI analysis unavailable',
                    transform=ax.transAxes, ha='center')
            ax.set_axis_off()
    fig.suptitle(f'PPI vs no PPI at the SAME distance setting | {level}')
    fig.savefig(out/'ppi_matched.png', dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6), constrained_layout=True)
    if len(paired):
        # r0_shape near 1h (v1 runs), else the stage-one working setting r0_raw near 2h.
        chosen = paired.iloc[0]
        for family, target in [('r0_shape', 1.0), ('r0_raw', 2.0)]:
            candidates = paired[paired.family == family]
            if len(candidates):
                chosen = candidates.iloc[(candidates.r0_h-target).abs().argmin()]
                break
        data = changes[changes.ppi_case == chosen.ppi_case].copy()
        data = data.iloc[np.argsort(-np.abs(data.delta_T.to_numpy()))[:12]]
        y = np.arange(len(data))
        ax.barh(y-.19, data.baseline_T, height=.34, label='no PPI', color=PALETTE['series2'])
        ax.barh(y+.19, data.ppi_T, height=.34, label='PPI', color=PALETTE['series1'])
        ax.set_yticks(y, [f'{r.ligand} -> {r.receptor} ({r.type})' for r in data.itertuples()])
        ax.set_xlabel('Transport amount summed over spatial connections')
        r0_text = f'r0/h={chosen.r0_h:g}' if pd.notna(chosen.r0_h) else 'd_ref/d'
        ax.set_title(f'{chosen.family}, {r0_text}; illustrative setting, not optimum')
        ax.legend(frameon=False)
        style_axis(ax)
    else:
        ax.text(.5, .5, 'No PPI input', transform=ax.transAxes, ha='center')
        ax.set_axis_off()
    fig.savefig(out/'lr_pair_changes.png', dpi=180)
    plt.close(fig)


def analyze(results, inputs, output):
    results, inputs, output = Path(results), Path(inputs), Path(output)
    manifest = json.loads((results/'manifest.json').read_text())
    check(manifest.get('status') == 'complete', 'Incomplete experiment')
    for filename, expected in manifest['input_sha256'].items():
        check(Path(filename).name == filename, 'Invalid input filename')
        check(digest(inputs/filename) == expected, f'Input changed since experiment: {filename}')
    summary = pd.read_csv(results/'summary.csv')
    check(summary.case.is_unique, 'Duplicate experiment cases')
    check(not summary.duplicated(['family','r0','score_mode']).any(), 'Duplicate conditions')
    pairs = pd.read_csv(results/'used_lr.csv', dtype=str)
    check(not pairs.duplicated(['ligand','receptor']).any(), 'Duplicate LR pairs')
    spatial = numeric_table(inputs/'spatial.csv')
    h = float(manifest['h'])
    check(np.isfinite(h) and h > 0, 'Invalid reference spacing')
    output.mkdir(parents=True, exist_ok=False)
    type_rows, lr_rows, paired_rows, changes, duplicate_rows = [], [], [], [], []
    control_rows, shared_rows = [], []
    artifact_hashes = {f: digest(results/f) for f in ['summary.csv','used_lr.csv','manifest.json']}
    canonical_ids = None

    def read_plan(case):
        nonlocal canonical_ids
        check(Path(case).name == case, 'Invalid case filename')
        filename = case+'.npz'
        with np.load(results/filename, allow_pickle=False) as z:
            ids = tuple(tuple(z[k].astype(str)) for k in ['spots','ligands','receptors'])
            if canonical_ids is None:
                canonical_ids = ids
            check(ids == canonical_ids, 'Spot/entity order changed between cases')
            spots, ligands, receptors = ids
            check(len(set(spots)) == len(spots) and set(spots) == set(spatial.index),
                  'Spatial spot IDs differ from saved plan')
            T, C, a, b = z['T'], z['C'], z['a'], z['b']
            check(T.shape == (len(a),len(b)), 'Budget/plan mismatch')
            check((T.sum(axis=1) <= a + 1e-7*max(1.,a.sum())).all()
                  and (T.sum(axis=0) <= b + 1e-7*max(1.,b.sum())).all(), 'Budget violation')
        artifact_hashes[filename] = digest(results/filename)
        return T, C, a, b, spots, ligands, receptors

    for (family, r0), group in summary.groupby(['family','r0'], dropna=False, sort=False):
        reference = group[group.score_mode == 'no_ppi']
        check(len(reference) == 1, 'Need exactly one no_ppi reference per distance setting')
        ref_case = reference.iloc[0]['case']
        base, _, base_a, base_b, _, _, _ = read_plan(ref_case)
        fingerprints = {}
        pair_maps = {}
        for row in group.itertuples(index=False):
            T, C, a, b, spots, ligands, receptors = read_plan(row.case)
            check(np.array_equal(a,base_a) and np.array_equal(b,base_b), 'Budgets changed between matched cases')
            check(np.isclose(T.sum(), row.total_T, atol=1e-8), 'Summary/plan total mismatch')
            xy = spatial.loc[list(spots),['x','y']].to_numpy(float)
            D = cdist(xy, xy)
            typed, per_lr = type_and_pair_metrics(T,C,spots,ligands,receptors,pairs,D)
            meta = dict(case=row.case, family=family, r0=r0, r0_h=r0/h, score_mode=row.score_mode)
            type_rows.extend([{**meta,**v} for v in typed])
            lr_rows.extend([{**meta,**v} for v in per_lr])
            pair_maps[row.case] = {(v['ligand'],v['receptor'],v['type']):v['total_T'] for v in per_lr}
            fingerprint = hashlib.sha256(np.ascontiguousarray(T).tobytes()).hexdigest()
            if fingerprint in fingerprints:
                duplicate_rows.append(dict(case=row.case, identical_plan_to=fingerprints[fingerprint],
                                           family=family, r0=r0))
            else:
                fingerprints[fingerprint] = row.case
            if row.score_mode == 'ppi':
                paired_rows.append(dict(family=family,r0=r0,r0_h=r0/h,
                    baseline_case=ref_case,ppi_case=row.case,**paired_metrics(base,T)))
            n = len(spots)
            li, ri = {v:i for i,v in enumerate(ligands)}, {v:i for i,v in enumerate(receptors)}
            blocks = [(slice(li[p.ligand]*n,(li[p.ligand]+1)*n), slice(ri[p.receptor]*n,(ri[p.receptor]+1)*n))
                      for p in pairs.itertuples(index=False)]
            if row.score_mode != 'no_ppi':
                control_rows.append(dict(family=family, r0=r0, r0_h=r0/h, baseline_case=ref_case,
                    case=row.case, score_mode=row.score_mode,
                    permuted_type=getattr(row, 'permuted_type', None),
                    **paired_metrics(base,T), **decompose(base,T,blocks)))
            masses = {(p.ligand,p.receptor): float(T[r,c].sum())
                      for p,(r,c) in zip(pairs.itertuples(index=False), blocks)}
            for role, column in [('ligand',0),('receptor',1)]:
                for entity in {k[column] for k in masses}:
                    partners = {k:v for k,v in masses.items() if k[column] == entity}
                    if len(partners) < 2:
                        continue
                    total = sum(partners.values())
                    for (ligand,receptor), mass in partners.items():
                        shared_rows.append(dict(case=row.case, family=family, r0=r0, r0_h=r0/h,
                            score_mode=row.score_mode, role=role, entity=entity, ligand=ligand,
                            receptor=receptor, mass=mass, share=mass/total if total else None))
        for row in group[group.score_mode == 'ppi'].itertuples(index=False):
            for (ligand,receptor,kind), value in pair_maps[row.case].items():
                previous = pair_maps[ref_case][ligand,receptor,kind]
                changes.append(dict(family=family,r0=r0,r0_h=r0/h,baseline_case=ref_case,
                    ppi_case=row.case,ligand=ligand,receptor=receptor,type=kind,
                    baseline_T=previous,ppi_T=value,delta_T=value-previous))

    types, per_lr = pd.DataFrame(type_rows), pd.DataFrame(lr_rows)
    paired = pd.DataFrame(paired_rows, columns=['family','r0','r0_h','baseline_case','ppi_case',
        'baseline_total_T','ppi_total_T','total_change_percent','paired_relative_L1','allocation_TV'])
    change_table = pd.DataFrame(changes, columns=['family','r0','r0_h','baseline_case','ppi_case',
        'ligand','receptor','type','baseline_T','ppi_T','delta_T'])
    controls = pd.DataFrame(control_rows, columns=['family','r0','r0_h','baseline_case','case','score_mode',
        'permuted_type','baseline_total_T','ppi_total_T','total_change_percent','paired_relative_L1',
        'allocation_TV','between_pair_TV','within_pair_TV'])
    shared = pd.DataFrame(shared_rows, columns=['case','family','r0','r0_h','score_mode','role','entity',
        'ligand','receptor','mass','share'])
    for name, frame in [('by_type.csv',types),('by_lr.csv',per_lr),
                        ('ppi_matched.csv',paired),('lr_pair_changes.csv',change_table),
                        ('control_matched.csv',controls),('shared_entity_split.csv',shared)]:
        frame.to_csv(output/name,index=False)
    pd.DataFrame(duplicate_rows, columns=['case','identical_plan_to','family','r0']).to_csv(
        output/'identical_plans.csv',index=False)
    cfg = manifest['config']
    plot_results(types, paired, change_table, output, cfg['coordinate_unit'], cfg['expression_level'])
    (output/'analysis_manifest.json').write_text(json.dumps(dict(
        status='complete', cases=len(summary), matched_ppi_comparisons=len(paired),
        reference='no_ppi with identical family and r0',
        expression_level=cfg['expression_level'], score_provenance=cfg.get('score_provenance'),
        artifact_sha256=artifact_hashes,
        notes=['Identical plans do not by themselves prove identical input scores.',
               'No r0 is selected as biologically optimal.',
               'Per-type unused budgets are not reported: entities can participate in both types.',
               'allocation_TV is missing when either total transport is zero.',
               'control_matched.csv: every score mode vs no_ppi; between_pair_TV = change of LR-pair '
               'shares, within_pair_TV = share-weighted change of each pair\'s own spatial plan.',
               'shared_entity_split.csv: how a ligand/receptor used by >=2 pairs splits its transport.']), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--results', default='results_demo')
    parser.add_argument('--input', default='demo_input')
    parser.add_argument('--output', default='analysis_demo')
    args = parser.parse_args()
    analyze(args.results, args.input, args.output)
