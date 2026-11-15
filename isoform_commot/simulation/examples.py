"""Reproduce COMMOT's PDE examples (Extended Data Fig. 1a, b) and check against the archived simulations.

  python -m simulation.examples --archive data/zenodo_7272562/pde_model_comparison.zip --output reports/pde_reimplementation

archive  For archived benchmark simulations (Zenodo 10.5281/zenodo.7272562), rerun the model from the archived
         inputs (ligand production, receptor initial field, rates, 50 x 50 grid) to t = 50 with
         (i) scipy RK45 at solve_ivp's default tolerances (rtol 1e-3, atol 1e-6; py-pde's ScipySolver default)
         and (ii) tight tolerances, and compare the complexes with the archived result.h5 at every saved time.
         Also compare t = 50 with the steady state.
1d       Two ligands on either side of one receptor (Extended Data Fig. 1a): unbiased; biased locations;
         higher diffusivity for ligand 2. Steady state, archived rates (D 200, a 1, b 0.001, c 0.2).
2d       Two crescent ligand sources, three receptors, pairs L1-R1, L1-R2, L2-R2, L2-R3 (Extended Data Fig. 1b).
No optimal transport is run here.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import io
import json
from pathlib import Path
import pickle
import time
import zipfile

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from research import check, launch_info
from simulation.pde import Grid, Model, System, bump, crescent, gaussian2d

RATES = dict(D=200.0, a=1.0, b=1e-3, c=0.2)            # identical in every archived simulation
ALLOWED = {('numpy.core.multiarray', '_reconstruct'), ('numpy', 'ndarray'), ('numpy', 'dtype')}


class _NumpyOnly(pickle.Unpickler):
    """meta_data.pkl holds numpy arrays and numbers only; refuse any other global."""
    def find_class(self, module, name):
        if (module, name) not in ALLOWED:
            raise pickle.UnpicklingError(f'blocked {module}.{name}')
        return getattr(importlib.import_module(module), name)


def load_archived(zf, case, sim):
    import h5py
    base = f'pde_model_comparison/data/case_{case}/sim_{sim}/pde/'
    meta = _NumpyOnly(io.BytesIO(zf.read(base + 'meta_data.pkl'))).load()
    a = np.asarray(meta['binding_rate'], float)
    nl, nr = a.shape
    npy = lambda name: np.load(io.BytesIO(zf.read(base + name)), allow_pickle=False)
    with h5py.File(io.BytesIO(zf.read(base + 'result.h5')), 'r') as h:
        data, times = h['data'][:], h['times'][:]
    labels = zf.read(base + 'result_labels.txt').decode().split()
    grid_xy = npy('grid.npy')
    model = Model(D=meta['ligand_diffusivity'], a=a, b=meta['dissociation_rate'], c=meta['degradation_rate_ligand'],
                  e=meta['degradation_rate_complex'], g=meta['degradation_rate_receptor'])
    lo, hi = meta['xlimit']
    grid = Grid(grid_xy.shape[0], float(lo), float(hi), 2)
    check(np.allclose(grid.coords(), grid_xy), 'Archived grid differs from a cell-centred grid on xlimit')
    system = System(model, grid, [npy(f'ligand_prod_{i+1}.npy') for i in range(nl)],
                    [npy(f'receptor_init_{j+1}.npy') for j in range(nr)])
    return system, data, times, labels


def archive_conservation(zf):
    """Receptor + complex totals of every archived simulation (conserved exactly by Eq. 19 with e = g = 0),
    with the solver start date and step count recorded by py-pde."""
    import h5py
    rows = []
    for case in range(1, 11):
        for sim in range(1, 11):
            base = f'pde_model_comparison/data/case_{case}/sim_{sim}/pde/'
            with h5py.File(io.BytesIO(zf.read(base + 'result.h5')), 'r') as h:
                data = h['data'][:]
                ctrl, solver = json.loads(h.attrs['controller']), json.loads(h.attrs['solver'])
            labels = zf.read(base + 'result_labels.txt').decode().split()
            total = sum(data[:, k].sum(axis=(1, 2)) for k, lab in enumerate(labels) if lab.startswith(('Receptor', 'Complex')))
            rows.append(dict(case=case, sim=sim, max_rel_drift=float(np.abs(total - total[0]).max() / total[0]),
                             solver_start=ctrl['solver_start'], steps=solver.get('steps'), t_final=ctrl.get('t_final'),
                             stop_reason=ctrl.get('stop_reason')))
    return pd.DataFrame(rows)


def complexes_of(system, fields_list):
    return np.stack([f[2] for f in fields_list])                  # times x pairs x grid


def archive_check(zf, cases):
    rows = []
    for case, sim in cases:
        s, data, times, labels = load_archived(zf, case, sim)
        archived = np.stack([data[:, labels.index(f'Complex_{i+1}_{j+1}')] for i, j in s.m.pairs], axis=1)
        receptors = np.stack([data[:, labels.index(f'Receptor_{j+1}')] for j in range(s.m.nr)], axis=1)
        drift = float(np.abs((receptors.sum(axis=(1, 2, 3)) + archived.sum(axis=(1, 2, 3)))
                             - (receptors[0].sum() + archived[0].sum())).max() / receptors[0].sum())
        start = time.time()
        loose, _ = s.integrate(times, rtol=1e-3, atol=1e-6)
        mid = time.time()
        tight, _ = s.integrate(times, rtol=1e-8, atol=1e-10)
        end = time.time()
        steady = s.steady_state()[2]
        for label, run in [('rk45_default_tol', loose), ('rk45_tight_tol', tight)]:
            ours = complexes_of(s, run)
            err = np.linalg.norm((ours - archived).reshape(len(times), -1), axis=1) / \
                np.maximum(np.linalg.norm(archived.reshape(len(times), -1), axis=1), 1e-300)
            final = ours[-1].ravel(), archived[-1].ravel()
            rows.append(dict(case=case, sim=sim, pairs=len(s.m.pairs), run=label,
                             rel_l2_t50=float(err[-1]), max_rel_l2_over_times=float(err[1:].max()),
                             spearman_t50=float(spearmanr(*final).statistic),
                             archived_receptor_complex_drift=drift,
                             our_receptor_complex_drift=float(np.abs(
                                 np.stack([f[1].sum() + f[2].sum() for f in run]) - s.r0.sum()).max() / s.r0.sum()),
                             rel_l2_t50_vs_steady=float(np.linalg.norm(ours[-1] - steady) / np.linalg.norm(steady)),
                             seconds=round((mid - start) if label.startswith('rk45_default') else (end - mid), 1)))
    return pd.DataFrame(rows)


def example_1d():
    grid = Grid(400, 0.0, 400.0, 1)
    x = grid.coords()
    setups = {'unbiased': dict(l1=300, l2=100, D=(200, 200)), 'biased_locations': dict(l1=270, l2=100, D=(200, 200)),
              'higher_D_ligand2': dict(l1=300, l2=100, D=(200, 800))}
    out = {}
    for name, s_ in setups.items():
        model = Model(D=s_['D'], a=[[RATES['a']], [RATES['a']]], b=[[RATES['b']], [RATES['b']]], c=[RATES['c']]*2)
        system = System(model, grid, [bump(x, s_['l1'], 50), bump(x, s_['l2'], 50)], [bump(x, 200, 50)])
        L, R, C = system.steady_state()
        _, _, Ct = system.integrate(np.array([0.0, 50.0]))[0][1]   # fields at t = 50
        out[name] = dict(system=system, setup=s_, L=L, R=R, C=C, C50=Ct)
    return x, out


def example_2d():
    grid = Grid(50, -200.0, 200.0, 2)
    xy = grid.coords()
    L1 = crescent(xy, (0, 0), 150, 18, -110, 110)                 # right-hand crescent
    L2 = crescent(xy, (0, 0), 150, 18, 120, 290)                  # left-hand crescent
    R1 = gaussian2d(xy, (-5, 30), (30, 16), 45)
    R2 = gaussian2d(xy, (15, -25), (28, 18), 100)
    R3 = gaussian2d(xy, (25, -10), (55, 14), 15)
    a = np.array([[1.0, 1.0, 0.0], [0.0, 1.0, 1.0]]) * RATES['a']
    model = Model(D=[RATES['D']]*2, a=a, b=np.where(a > 0, RATES['b'], 0.0), c=[RATES['c']]*2)
    system = System(model, grid, [L1, L2], [R1, R2, R3])
    L, R, C = system.steady_state()
    (_, _, C50) = system.integrate(np.array([0.0, 50.0]))[0][1]
    return xy, system, L, R, C, C50


def plot_1d(x, out, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from analyze import PALETTE, style_axis
    titles = {'unbiased': 'Unbiased diffusivity and locations', 'biased_locations': 'Biased locations',
              'higher_D_ligand2': 'Higher diffusivity for ligand 2'}
    fig, axes = plt.subplots(2, 3, figsize=(11, 4.6), sharex=True, sharey=True, constrained_layout=True)
    for col, (name, r) in enumerate(out.items()):
        s = r['system']
        for i in range(2):
            ax = axes[i, col]
            prod = s.p[i] / s.p[i].max()
            ax.fill_between(x, prod, color=PALETTE['series1'], alpha=.25, lw=0, label='ligand production (input)')
            ax.fill_between(x, s.r0[0], color=PALETTE['series2'], alpha=.25, lw=0, label='receptor (initial)')
            ax.fill_between(x, r['C'][i], color=PALETTE['series2'], alpha=.9, lw=0, label='occupied receptor, steady state')
            ax.plot(x, r['C50'][i], color=PALETTE['text'], lw=1, ls='--', label='occupied receptor, t = 50')
            ax.text(8, .85, f'L{i+1}', fontsize=9, color=PALETTE['text'])
            if i == 0:
                ax.set_title(titles[name], fontsize=9, loc='left')
            ax.set_ylim(0, 1.05)
            style_axis(ax)
    axes[1, 0].set_xlabel('position'); axes[1, 1].set_xlabel('position'); axes[1, 2].set_xlabel('position')
    axes[0, 0].legend(fontsize=7, frameon=False, loc='upper center')
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_2d(xy, system, C, C50, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ext = [system.grid.lo, system.grid.hi, system.grid.lo, system.grid.hi]
    names = [f'L{i+1}-R{j+1}' for i, j in system.m.pairs]
    fig, axes = plt.subplots(3, 5, figsize=(13, 8), constrained_layout=True)

    def show(ax, f, cmap, title, vmax):
        im = ax.imshow(f.T, origin='lower', extent=ext, cmap=cmap, vmin=0, vmax=vmax)
        ax.set_title(title, fontsize=9)
        return im
    inputs = [(system.p[i], 'Blues', f'L{i+1} production') for i in range(2)] + \
             [(system.r0[j], 'Reds', f'R{j+1} initial') for j in range(3)]
    for ax, (f, cmap, title) in zip(axes[0], inputs):
        show(ax, f, cmap, title, f.max())
    for row, fields, label in [(1, C, 'steady state'), (2, C50, 't = 50')]:
        vmax = fields.max()                                     # common scale within a row
        for k, name in enumerate(names):
            im = show(axes[row, k], fields[k], 'Reds', f'{name}, {label}', vmax)
        axes[row, 4].axis('off')
        fig.colorbar(im, ax=axes[row, 4], location='left', shrink=.8, label=f'complex ({label}; common scale)')
    for ax in axes.ravel():
        ax.set_xticks([]); ax.set_yticks([])
    fig.savefig(path, dpi=160)
    plt.close(fig)


def run(args):
    launch = launch_info()
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    zpath = Path(args.archive)
    check(zpath.is_file(), 'archive zip missing')
    digest = hashlib.sha256(zpath.read_bytes()).hexdigest()
    check(digest == args.archive_sha256, f'archive SHA256 {digest} differs from the recorded value')
    zf = zipfile.ZipFile(zpath)
    cases = [tuple(map(int, c.split(':'))) for c in args.cases]
    table = archive_check(zf, cases)
    conservation = archive_conservation(zf)
    x, ex1 = example_1d()
    xy, sys2, L2, R2, C2, C50_2 = example_2d()
    out.mkdir(parents=True)
    table.to_csv(out/'archive_check.csv', index=False)
    conservation.to_csv(out/'archive_conservation.csv', index=False)
    rows = []
    for name, r in ex1.items():
        for i in range(2):
            rows.append(dict(example=name, ligand=f'L{i+1}', D=r['setup']['D'][i], ligand_centre=r['setup'][f'l{i+1}'],
                             occupied_share=float(r['C'][i].sum() / r['system'].r0.sum()),
                             occupied_share_t50=float(r['C50'][i].sum() / r['system'].r0.sum()),
                             peak_position=float(x[np.argmax(r['C'][i])])))
    pd.DataFrame(rows).to_csv(out/'example_1d.csv', index=False)
    names = [f'L{i+1}-R{j+1}' for i, j in sys2.m.pairs]
    pd.DataFrame(dict(pair=names, total_steady=C2.reshape(len(names), -1).sum(1),
                      total_t50=C50_2.reshape(len(names), -1).sum(1),
                      share_of_receptor=[C2[k].sum() / sys2.r0[j].sum() for k, (i, j) in enumerate(sys2.m.pairs)])).to_csv(
        out/'example_2d.csv', index=False)
    plot_1d(x, ex1, out/'example_1d.png')
    plot_2d(xy, sys2, C2, C50_2, out/'example_2d.png')
    (out/'run.json').write_text(json.dumps(dict(status='complete', launch=launch,
        created=datetime.now(timezone.utc).isoformat(timespec='seconds'), archive=str(zpath), archive_sha256=digest,
        archive_cases=args.cases, rates=RATES, no_optimal_transport=True), indent=2))
    print(table.to_string(index=False))
    print(pd.DataFrame(rows).to_string(index=False))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--archive', default='data/zenodo_7272562/pde_model_comparison.zip')
    p.add_argument('--archive-sha256', default='ec8eb74c878f299b0b7c017a1ca758777c188c2ea1c672120e62b2de42ccc50f')
    p.add_argument('--cases', nargs='+', default=['1:1', '2:1', '6:1', '10:1'], help='case:sim pairs to rerun')
    p.add_argument('--output', required=True)
    return p


if __name__ == '__main__':
    run(parser().parse_args())
