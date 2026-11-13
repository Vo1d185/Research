"""Original COMMOT collective OT as a baseline on the SAME prepared input.

Calls the official cot_combine_sparse (COMMOT commit d117445) unchanged:
Euclidean distance cost, transport forbidden beyond dis_thr, cost normalized by
the largest eligible cost, and the COMMOT default mixture of four transport
setups (all-all, row, column, block; weights 0.25 each) with eps=0.1, rho=10.
Only the all-all setup respects every expression budget jointly; the mixture
can exceed a single entity's budget (recorded in summary.csv). --weights 1 0 0 0
runs the pure collective setup.
As in COMMOT, contact and diffusion pairs are NOT distinguished and same-spot
transport (d=0) is allowed. Output files follow research.py (summary.csv,
case_*.npz, used_lr.csv, manifest.json), so evaluate_targets.py and
analyze.py-style readers can use them directly.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from fetch_solver import COMMIT, COT_SHA256, COT_URL, PACKAGE, SHA256, URL, VENDOR
from research import check, numeric_table


def official_cot():
    for path, sha in [(PACKAGE/'_unot.py', SHA256), (PACKAGE/'_cot.py', COT_SHA256)]:
        check(path.is_file(), 'Run python fetch_solver.py first')
        check(hashlib.sha256(path.read_bytes()).hexdigest() == sha, f'Checksum mismatch: {path}')
    import matplotlib
    matplotlib.use('Agg')  # the official _unot.py imports pyplot
    if str(VENDOR) not in sys.path:
        sys.path.insert(0, str(VENDOR))
    return importlib.import_module('commot_ot._cot').cot_combine_sparse


def run(folder, output, dis_thr_h=(3.0, 5.0, 10.0), eps=0.1, rho=10.0,
        weights=(0.25, 0.25, 0.25, 0.25), nitermax=10000):
    folder, output = Path(folder), Path(output)
    check(not output.exists(), f'Output exists: {output}')
    cfg = json.loads((folder/'config.json').read_text())
    check(cfg.get('expression_state') == 'processed', 'Use a prepared (processed) input directory')
    check(len(weights) == 4 and abs(sum(weights) - 1) < 1e-9 and min(weights) >= 0, 'Invalid weights')
    check(len(dis_thr_h) > 0 and min(dis_thr_h) > 0, 'Invalid dis_thr_h')
    x, spatial = numeric_table(folder/'expression.csv'), numeric_table(folder/'spatial.csv')
    check(set(x.index) == set(spatial.index), 'Spot IDs must match')
    xy = spatial.loc[x.index, ['x', 'y']].to_numpy(float)
    lr = pd.read_csv(folder/'lr.csv', dtype=str)[['ligand', 'receptor', 'type']].drop_duplicates()
    positive = set(x.columns[x.sum() > 0])
    lr = lr[lr.ligand.isin(positive) & lr.receptor.isin(positive)].reset_index(drop=True)
    check(len(lr) > 0, 'No expressed LR pairs')
    ligands, receptors = sorted(lr.ligand.unique()), sorted(lr.receptor.unique())
    n = len(x)
    check(len(ligands)*len(receptors)*n*n <= 4000000, 'Dense output too large for this prototype')
    D = cdist(xy, xy)
    nearest = D.copy(); np.fill_diagonal(nearest, np.inf)
    h = float(np.median(nearest.min(axis=1)))
    S, R = x[ligands].to_numpy(float), x[receptors].to_numpy(float)
    A = np.full((len(ligands), len(receptors)), np.inf)
    for row in lr.itertuples(index=False):
        A[ligands.index(row.ligand), receptors.index(row.receptor)] = 1.0
    a, b = S.T.reshape(-1), R.T.reshape(-1)
    cot = official_cot()
    output.mkdir(parents=True)
    lr.to_csv(output/'used_lr.csv', index=False)
    records = []
    for k, multiple in enumerate(sorted(set(dis_thr_h))):
        dis_thr = multiple*h
        P = cot(S, R, A, D, dis_thr*np.ones_like(A), eps_p=eps, rho=rho, weights=tuple(weights),
                nitermax=nitermax)
        T = np.zeros((len(ligands)*n, len(receptors)*n))
        C = np.full_like(T, np.inf)
        for (i, j), plan in P.items():
            T[i*n:(i+1)*n, j*n:(j+1)*n] = plan.toarray()
            block = C[i*n:(i+1)*n, j*n:(j+1)*n]
            block[D <= dis_thr] = D[D <= dis_thr]  # eligibility record (unnormalized distance)
        check(np.isfinite(T).all() and (T >= -1e-12).all(), 'Invalid COMMOT transport')
        T = np.maximum(T, 0)
        check((T[~np.isfinite(C)] == 0).all(), 'Transport beyond dis_thr')
        # The default mixture is NOT budget-feasible by design: the row, column and
        # block setups solve sub-problems that each reuse the full budgets. Record it.
        excess_l, excess_r = float((T.sum(axis=1) - a).max()), float((T.sum(axis=0) - b).max())
        name = f'case_{k:03d}'
        np.savez_compressed(output/f'{name}.npz', C=C, T=T, a=a, b=b, spots=x.index.to_numpy(str),
            ligands=np.array(ligands), receptors=np.array(receptors), scores=np.ones(len(lr)),
            pair_ligand=lr.ligand.to_numpy(str), pair_receptor=lr.receptor.to_numpy(str),
            pair_type=lr.type.to_numpy(str))
        spot = T.reshape(len(ligands), n, len(receptors), n).sum(axis=(0, 2))
        total = float(T.sum())
        records.append(dict(case=name, family='commot_euc', r0=dis_thr, dis_thr_h=multiple,
            score_mode='no_ppi', permuted_type=None, total_T=total, unused_ligand=float(a.sum()-total),
            unused_receptor=float(b.sum()-total), max_ligand_budget_excess=excess_l,
            max_receptor_budget_excess=excess_r,
            mean_distance=float((spot*D).sum()/total) if total else None))
        print(f'COMMOT dis_thr={multiple:g}h total_T={total:.4g}', flush=True)
    pd.DataFrame(records).to_csv(output/'summary.csv', index=False)
    files = ['config.json', 'expression.csv', 'spatial.csv', 'lr.csv']
    (output/'manifest.json').write_text(json.dumps(dict(status='complete', method='COMMOT original',
        config=dict(coordinate_unit=cfg['coordinate_unit'], expression_level=cfg['expression_level'],
                    eps=eps, rho=rho, weights=list(weights), nitermax=nitermax,
                    dis_thr_h=sorted(set(dis_thr_h))), h=h,
        solver=dict(commit=COMMIT, files={URL: SHA256, COT_URL: COT_SHA256}, function='cot_combine_sparse'),
        budget_feasible=bool(weights[1] == weights[2] == weights[3] == 0),
        notes=['The default 4-setup mixture can exceed per-entity expression budgets by design '
               '(row/column/block setups reuse budgets); --weights 1 0 0 0 gives the pure collective setup.',
               'Contact/diffusion types not distinguished, as in COMMOT.',
               'Same-spot transport allowed (d=0), unlike research.py.',
               'C in case_*.npz records eligibility with raw distances; COMMOT normalizes internally.'],
        input_sha256={f: hashlib.sha256((folder/f).read_bytes()).hexdigest() for f in files}), indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True); p.add_argument('--output', required=True)
    p.add_argument('--dis-thr-h', type=float, nargs='+', default=[3.0, 5.0, 10.0],
                   help='distance cutoffs as multiples of h (default 3 5 10; 5h is the main setting)')
    p.add_argument('--eps', type=float, default=0.1); p.add_argument('--rho', type=float, default=10.0)
    p.add_argument('--weights', type=float, nargs=4, default=[0.25, 0.25, 0.25, 0.25])
    args = p.parse_args()
    run(args.input, args.output, args.dis_thr_h, args.eps, args.rho, args.weights)
