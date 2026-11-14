"""Original COMMOT collective OT as a baseline on the SAME prepared input.

Calls the official cot_combine_sparse (COMMOT commit d117445) unchanged:
Euclidean distance cost, transport forbidden beyond dis_thr, cost normalized by
the largest eligible cost, and the COMMOT default mixture of four transport
setups (all-all, row, column, block; weights 0.25 each) with eps=0.1, rho=10.
Only the all-all setup respects every expression budget jointly; the mixture
can exceed a single entity's budget (recorded in summary.csv). --weights 1 0 0 0
runs the pure collective setup: it calls COMMOT's cot_sparse directly (the collective
setup inside cot_combine_sparse; the other three setups would get weight 0) and applies the
same mass-conservation check as research.solve_sparse, using COMMOT's internally normalized
cost d / cost_scale (cost_scale = largest eligible distance).
As in COMMOT, contact and diffusion pairs are NOT distinguished and same-spot
transport (d=0) is allowed. Output files follow research.py (summary.csv,
sparse COO case_*.npz via plans.py, used_lr.csv, manifest.json), so evaluate_targets.py and
analyze.py-style readers can use them directly.
"""
import argparse
import hashlib
import importlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist

from fetch_solver import COMMIT, COT_SHA256, COT_URL, PACKAGE, SHA256, URL, VENDOR
import plans
from research import check, numeric_table, sparse_kkt_check, sparse_mass_check


def official_cot():
    for path, sha in [(PACKAGE/'_unot.py', SHA256), (PACKAGE/'_cot.py', COT_SHA256)]:
        check(path.is_file(), 'Run python fetch_solver.py first')
        check(hashlib.sha256(path.read_bytes()).hexdigest() == sha, f'Checksum mismatch: {path}')
    import matplotlib
    matplotlib.use('Agg')  # the official _unot.py imports pyplot
    if str(VENDOR) not in sys.path:
        sys.path.insert(0, str(VENDOR))
    module = importlib.import_module('commot_ot._cot')
    return module.cot_combine_sparse, module.cot_sparse


def run(folder, output, dis_thr_h=(3.0, 5.0, 10.0), eps=0.1, rho=10.0,
        weights=(0.25, 0.25, 0.25, 0.25), nitermax=10000, dis_thr_um=None, stopthr=1e-8):
    """dis_thr_um (micrometre inputs) overrides dis_thr_h. Case files are sparse COO (plans.py).
    nitermax/stopthr default to COMMOT's own defaults (1e4, 1e-8)."""
    folder, output = Path(folder), Path(output)
    check(not output.exists(), f'Output exists: {output}')
    cfg = json.loads((folder/'config.json').read_text())
    check(cfg.get('expression_state') == 'processed', 'Use a prepared (processed) input directory')
    weights = tuple(float(w) for w in weights)
    check(len(weights) == 4 and abs(sum(weights) - 1) < 1e-9 and min(weights) >= 0, 'Invalid weights')
    pure = weights == (1.0, 0.0, 0.0, 0.0)
    if dis_thr_um:
        check(cfg.get('coordinate_unit') == 'micrometer', '--dis-thr-um needs micrometre coordinates')
        check(min(dis_thr_um) > 0, 'Invalid dis_thr_um')
    else:
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
    D = cdist(xy, xy)
    nearest = D.copy(); np.fill_diagonal(nearest, np.inf)
    h = float(np.median(nearest.min(axis=1)))
    S, R = x[ligands].to_numpy(float), x[receptors].to_numpy(float)
    A = np.full((len(ligands), len(receptors)), np.inf)
    for row in lr.itertuples(index=False):
        A[ligands.index(row.ligand), receptors.index(row.receptor)] = 1.0
    a, b = S.T.reshape(-1), R.T.reshape(-1)
    combine, collective = official_cot()
    max_amount = max(a.sum(), b.sum())   # COMMOT cot_sparse: budgets / max(S.sum(), D.sum())
    cutoffs = [(float(c), float(c)/h) for c in sorted(set(dis_thr_um))] if dis_thr_um else \
              [(m*h, float(m)) for m in sorted(set(dis_thr_h))]
    output.mkdir(parents=True)
    lr.to_csv(output/'used_lr.csv', index=False)
    records = []
    for k, (dis_thr, multiple) in enumerate(cutoffs):
        cut = dis_thr*np.ones_like(A)
        start = time.time()
        if pure:
            P = collective(S, R, A, D, cut, eps_p=eps, rho=rho, nitermax=nitermax, stopthr=stopthr)
        else:
            P = combine(S, R, A, D, cut, eps_p=eps, rho=rho, weights=weights, nitermax=nitermax, stopthr=stopthr)
        seconds = time.time() - start
        rows, cols, Ts, Cs = [], [], [], []
        for (i, j), plan in P.items():
            plan = sparse.coo_matrix(plan)
            plan.sum_duplicates()
            rows.append(i*n + plan.row); cols.append(j*n + plan.col)
            Ts.append(plan.data); Cs.append(D[plan.row, plan.col])
        row, col = np.concatenate(rows).astype(np.int64), np.concatenate(cols).astype(np.int64)
        T, C = np.concatenate(Ts), np.concatenate(Cs)  # C: eligibility record (raw distance)
        check(np.isfinite(T).all() and (T >= -1e-12).all(), 'Invalid COMMOT transport')
        T = np.maximum(T, 0)
        check((C <= dis_thr*(1 + 1e-12)).all(), 'Transport beyond dis_thr')
        check(len(np.unique(row*len(b) + col)) == len(row), 'Duplicate entries in COMMOT output')
        # The default mixture is NOT budget-feasible by design: the row, column and
        # block setups solve sub-problems that each reuse the full budgets. Record it.
        excess_l = float((np.bincount(row, weights=T, minlength=len(a)) - a).max())
        excess_r = float((np.bincount(col, weights=T, minlength=len(b)) - b).max())
        cost_scale = float(D[D <= dis_thr].max())   # as in COMMOT cot_sparse
        checked = {}
        if pure:   # same check as research.solve_sparse, on COMMOT's normalized problem
            try:
                args = (row, col, C/cost_scale, T/max_amount, a/max_amount, b/max_amount, eps, rho)
                info = dict(sparse_mass_check(*args), **sparse_kkt_check(*args))
                checked = dict(check_passed=True, **{f'check_{key}': v for key, v in info.items()})
            except ValueError as error:
                checked = dict(check_passed=False, check_error=str(error))
                print(f'WARNING: mass check failed at dis_thr={dis_thr:g}: {error}', flush=True)
        name = f'case_{k:03d}'
        plans.save_sparse(output/f'{name}.npz', row, col, C, T, (len(ligands)*n, len(receptors)*n), a, b,
            x.index.to_numpy(str), ligands, receptors, scores=np.ones(len(lr)),
            pair_ligand=lr.ligand.to_numpy(str), pair_receptor=lr.receptor.to_numpy(str),
            pair_type=lr.type.to_numpy(str))
        total = float(T.sum())
        records.append(dict(case=name, family='commot_euc', r0=dis_thr, dis_thr_h=multiple,
            dis_thr_um=dis_thr if dis_thr_um else None,
            score_mode='no_ppi', permuted_type=None, total_T=total, unused_ligand=float(a.sum()-total),
            unused_receptor=float(b.sum()-total), max_ligand_budget_excess=excess_l,
            max_receptor_budget_excess=excess_r, entries=len(T),
            mean_distance=float((T*C).sum()/total) if total else None, eps=eps, rho=rho,
            scheme='collective' if pure else 'mixture', weights=' '.join(f'{w:g}' for w in weights),
            nitermax=nitermax, stopthr=stopthr, cost_scale=cost_scale, frac_ligand=total/float(a.sum()),
            frac_receptor=total/float(b.sum()), solve_seconds=round(seconds, 1), **checked))
        print(f'COMMOT dis_thr={dis_thr:g} ({multiple:g}h) total_T={total:.4g} {seconds:.0f}s', flush=True)
    pd.DataFrame(records).to_csv(output/'summary.csv', index=False)
    files = ['config.json', 'expression.csv', 'spatial.csv', 'lr.csv']
    (output/'manifest.json').write_text(json.dumps(dict(status='complete', method='COMMOT original',
        config=dict(coordinate_unit=cfg['coordinate_unit'], expression_level=cfg['expression_level'],
                    eps=eps, rho=rho, weights=list(weights), nitermax=nitermax, stopthr=stopthr,
                    dis_thr_h=[c[1] for c in cutoffs], dis_thr_um=[c[0] for c in cutoffs] if dis_thr_um else None),
        h=h, format='sparse coo (plans.py)',
        solver=dict(commit=COMMIT, files={URL: SHA256, COT_URL: COT_SHA256},
                    function='cot_sparse' if pure else 'cot_combine_sparse'),
        budget_feasible=bool(weights[1] == weights[2] == weights[3] == 0),
        notes=['The default 4-setup mixture can exceed per-entity expression budgets by design '
               '(row/column/block setups reuse budgets); --weights 1 0 0 0 gives the pure collective setup.',
               'Contact/diffusion types not distinguished, as in COMMOT.',
               'Same-spot transport allowed (d=0), unlike stage-1 research.py.',
               'Stored entries are those COMMOT returns (within dis_thr, nonzero sender and receiver mass); '
               'C records the raw distance; COMMOT normalizes the cost internally.'],
        input_sha256={f: hashlib.sha256((folder/f).read_bytes()).hexdigest() for f in files}), indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input', required=True); p.add_argument('--output', required=True)
    p.add_argument('--dis-thr-h', type=float, nargs='+', default=[3.0, 5.0, 10.0],
                   help='distance cutoffs as multiples of h (default 3 5 10; 5h is the main setting)')
    p.add_argument('--eps', type=float, default=0.1); p.add_argument('--rho', type=float, default=10.0)
    p.add_argument('--weights', type=float, nargs=4, default=[0.25, 0.25, 0.25, 0.25])
    p.add_argument('--nitermax', type=int, default=10000, help='COMMOT default 1e4')
    p.add_argument('--stopthr', type=float, default=1e-8, help='COMMOT default 1e-8')
    p.add_argument('--dis-thr-um', type=float, nargs='+',
                   help='cutoffs in micrometres (micrometre inputs, e.g. stage 2); overrides --dis-thr-h')
    args = p.parse_args()
    run(args.input, args.output, args.dis_thr_h, args.eps, args.rho, args.weights, args.nitermax,
        dis_thr_um=args.dis_thr_um, stopthr=args.stopthr)
