"""GPU (torch) port of COMMOT's unot_sinkhorn_l1_sparse, for reports/stage2_inject_01 (amendment_01.md).

The iteration is a line-by-line transcription of vendor/_unot.py::unot_sinkhorn_l1_sparse (COMMOT commit
d117445b): same updates, same relative-change stopping rule checked every 10 iterations, same returned
iteration count (loop count; > nitermax means the cap was reached), float64 throughout. solve_sparse_gpu
mirrors research.solve_sparse (active rows/columns, common scale) and runs the SAME mass and KKT checks on
the result (research.sparse_mass_check, research.sparse_kkt_check), so convergence is judged identically.

  equivalence  build one stage2_inject cost, run the official CPU solver and this port for a fixed number
               of iterations (stopthr = 0) on the identical normalised problem, report the largest
               relative difference of the transport entries. run_inject.sh refuses GPU runs above 1e-9.

Environment: any env with torch (CUDA build for --device cuda), numpy, scipy, pandas, plus this repository.
"""
import argparse
import json
from pathlib import Path

import numpy as np

EQUIVALENCE_TOL = 1e-9


def torch_device(device):
    """'cuda' / 'cuda:N' run the port on a GPU; 'torch-cpu' runs the same port on CPU (local tests only)."""
    return 'cpu' if device == 'torch-cpu' else device


def unot_sinkhorn_l1_sparse_torch(a, b, row, col, c, eps, m, nitermax=10000, stopthr=1e-8, device='cuda'):
    """Same computation as the official solver; returns (transport on the given entries, iteration count)."""
    import torch
    device = torch_device(device)
    torch.use_deterministic_algorithms(True, warn_only=True)
    dt = torch.float64
    a_ = torch.as_tensor(np.asarray(a, float), dtype=dt, device=device)
    b_ = torch.as_tensor(np.asarray(b, float), dtype=dt, device=device)
    r = torch.as_tensor(np.asarray(row, np.int64), device=device)
    s = torch.as_tensor(np.asarray(col, np.int64), device=device)
    C = torch.as_tensor(np.asarray(c, float), dtype=dt, device=device)
    f = torch.zeros_like(a_); g = torch.zeros_like(b_)
    loga, logb = torch.log(a_), torch.log(b_)
    niter, err = 0, 100.0
    while niter <= nitermax and err > stopthr:
        fprev, gprev = f, g
        K = torch.exp((-C + f[r] + g[s]) / eps)
        rows = torch.zeros_like(a_).index_add_(0, r, K)
        f = eps * loga - eps * torch.log(rows + torch.exp((-m + f) / eps)) + f
        K = torch.exp((-C + f[r] + g[s]) / eps)
        cols = torch.zeros_like(b_).index_add_(0, s, K)
        g = eps * logb - eps * torch.log(cols + torch.exp((-m + g) / eps)) + g
        if niter % 10 == 0:
            err_f = (f - fprev).abs().max() / torch.maximum(torch.maximum(f.abs().max(), fprev.abs().max()),
                                                            torch.ones((), dtype=dt, device=device))
            err_g = (g - gprev).abs().max() / torch.maximum(torch.maximum(g.abs().max(), gprev.abs().max()),
                                                            torch.ones((), dtype=dt, device=device))
            err = float(0.5 * (err_f + err_g))
        niter += 1
    K = torch.exp((-C + f[r] + g[s]) / eps)
    return K.cpu().numpy(), niter


def solve_sparse_gpu(row, col, c, a, b, eps=0.2, rho=1.0, nitermax=50000, stopthr=1e-11, device='cuda'):
    """research.solve_sparse with the solver call replaced by the torch port (same checks, same outputs)."""
    from research import check, sparse_mass_check, sparse_kkt_check
    a, b = np.asarray(a, float), np.asarray(b, float)
    row, col, c = np.asarray(row, np.int64), np.asarray(col, np.int64), np.asarray(c, float)
    check(a.ndim == b.ndim == 1 and row.shape == col.shape == c.shape, 'Wrong dimensions')
    check(np.isfinite(a).all() and np.isfinite(b).all() and (a >= 0).all() and (b >= 0).all(),
          'Invalid expression budgets')
    check(np.isfinite([eps, rho]).all() and min(eps, rho) > 0, 'Invalid OT parameters')
    check(len(c) == 0 or (np.isfinite(c).all() and (c >= 0).all()), 'Invalid costs')
    check(len(np.unique(row * len(b) + col)) == len(row), 'Duplicate entries')
    T = np.zeros(len(c))
    scale = max(float(a.sum()), float(b.sum()))
    keep = (a[row] > 0) & (b[col] > 0) if len(c) else np.zeros(0, bool)
    if not scale or not keep.any():
        return T, dict(entries=0)
    ia, r = np.unique(row[keep], return_inverse=True)
    ib, s = np.unique(col[keep], return_inverse=True)
    aa, bb = a[ia] / scale, b[ib] / scale
    K, niter = unot_sinkhorn_l1_sparse_torch(aa, bb, r, s, c[keep], eps, rho, nitermax, stopthr, device)
    check(np.isfinite(K).all(), 'Nonfinite transport')
    info = sparse_mass_check(r, s, c[keep], K, aa, bb, eps, rho)
    info.update(sparse_kkt_check(r, s, c[keep], K, aa, bb, eps, rho))
    T[keep] = K * scale
    info.update(entries=int(keep.sum()), active_rows=len(ia), active_cols=len(ib), niter=niter,
                nitermax=int(nitermax), cap_reached=bool(niter > nitermax), backend=f'torch:{device}')
    return T, info


def normalised_problem(row, col, c, a, b):
    """The exact (aa, bb, r, s, c) that research.solve_sparse hands to the solver."""
    keep = (a[row] > 0) & (b[col] > 0)
    scale = max(float(a.sum()), float(b.sum()))
    ia, r = np.unique(row[keep], return_inverse=True)
    ib, s = np.unique(col[keep], return_inverse=True)
    return a[ia] / scale, b[ib] / scale, r, s, c[keep]


def compare(aa, bb, r, s, c, eps, rho, iters, device):
    from scipy import sparse
    from research import backend_sparse
    import time
    C = sparse.coo_matrix((c, (r, s)), shape=(len(aa), len(bb)))
    t0 = time.perf_counter()
    K_cpu = backend_sparse()(aa, bb, C, eps, rho, nitermax=iters, stopthr=0.0)
    t1 = time.perf_counter()
    K_gpu, n = unot_sinkhorn_l1_sparse_torch(aa, bb, r, s, c, eps, rho, iters, 0.0, device)   # returns on host: synced
    t2 = time.perf_counter()
    check_order = np.array_equal(K_cpu.row, r) and np.array_equal(K_cpu.col, s)
    scale = float(np.abs(K_cpu.data).max()) or 1.0
    return dict(entries=int(len(c)), iterations=int(n), order_preserved=bool(check_order),
                max_abs_diff=float(np.abs(K_cpu.data - K_gpu).max()),
                max_rel_diff=float(np.abs(K_cpu.data - K_gpu).max() / scale),
                cpu_seconds=round(t1 - t0, 3), port_seconds=round(t2 - t1, 3),   # timing only; not a criterion
                speedup=round((t1 - t0) / max(t2 - t1, 1e-9), 2))


def equivalence(args):
    import pandas as pd
    import stage2
    import stage2_inject as si
    prep, cfg, x, spatial, lr, ligands, receptors, a, b, edges = si.load(args.input)
    row, col, c0, d, p = stage2.build_sparse_cost(lr, ligands, receptors, a, b, len(x), edges, args.r0_um)
    if args.variant == 'none':
        c = c0
    else:
        sc = pd.read_csv(Path(args.scores)/'scores.csv')
        c = c0 + args.beta * (-np.log(si.effective_score(lr, sc[args.variant], args.norm)))[p]
    res = compare(*normalised_problem(row, col, c, a, b), args.eps, args.rho, args.iters, args.device)
    res.update(variant=args.variant, norm=args.norm, beta=args.beta, eps=args.eps, rho=args.rho,
               device=args.device, tolerance=EQUIVALENCE_TOL, passed=bool(res['order_preserved'] and
                                                                          res['max_rel_diff'] <= EQUIVALENCE_TOL))
    out = Path(args.output)
    if out.exists():
        raise FileExistsError(f'{out} exists')
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2))
    print(json.dumps(res, indent=2))
    if not res['passed']:
        raise SystemExit('GPU port differs from the official solver')


if __name__ == '__main__':
    pa = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = pa.add_subparsers(dest='command', required=True)
    e = sub.add_parser('equivalence')
    for k in ('input', 'scores', 'output'):
        e.add_argument('--' + k, required=(k != 'scores'))
    e.add_argument('--variant', default='none'); e.add_argument('--norm', default='raw', choices=['raw', 'ligand'])
    e.add_argument('--beta', type=float, default=1.0); e.add_argument('--r0-um', type=float, default=50.0)
    e.add_argument('--eps', type=float, required=True); e.add_argument('--rho', type=float, required=True)
    e.add_argument('--iters', type=int, default=500); e.add_argument('--device', default='cuda')
    a_ = pa.parse_args()
    equivalence(a_)
