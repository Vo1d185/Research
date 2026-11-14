"""Experimental costs in the official collective L1-UNOT core (not full COMMOT)."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.spatial import Delaunay, QhullError
from fetch_solver import TARGET, SHA256, URL


def check(condition, message):
    if not condition:
        raise ValueError(message)


def _official_unot():
    check(TARGET.is_file(), 'Run python fetch_solver.py first')
    check(hashlib.sha256(TARGET.read_bytes()).hexdigest() == SHA256,
          'Official solver checksum mismatch')
    spec = importlib.util.spec_from_file_location('official_unot', TARGET)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def backend():
    return _official_unot().unot_sinkhorn_l1_dense


def backend_sparse():
    return _official_unot().unot_sinkhorn_l1_sparse


def sparse_mass_check(row, col, c, t, aa, bb, eps, rho, tol=1e-7):
    """Mass balance for the sparse solver, which returns T but not the duals.

    The official iteration ends with a g-update, so column constraints hold exactly:
    nu = b - T^T 1. At the optimum T_ij = mu_i nu_j exp((2 rho - C_ij)/eps), so mu is
    recovered from T and nu; then |T 1 + mu - a| < tol is the same row check as the
    dense solver's, and every stored entry must have that multiplicative form.
    All quantities are on the normalized scale (budgets divided by the common scale).

    nu is only resolvable where the unused receptor mass exceeds rounding error; with a
    large rho/eps it can fall far below double precision (e.g. 1e-88), where b - T^T 1 is
    noise. Rows are checked through resolvable columns only; the numbers of checked rows
    and columns are reported. Feasibility (no budget exceeded) is checked everywhere.
    """
    check(((2*rho - c.min()) / eps) < 700, 'exp((2 rho - C)/eps) overflows; use a smaller rho/eps')
    E = np.exp((2*rho - c) / eps)
    rowsum = np.bincount(row, weights=t, minlength=len(aa))
    colsum = np.bincount(col, weights=t, minlength=len(bb))
    nu = bb - colsum
    check(nu.min() > -tol and (aa - rowsum).min() > -tol, 'Budget exceeded: solver has not converged')
    ok_col = nu > 1e-9 * bb
    use = ok_col[col]
    num = np.bincount(row, weights=np.where(use, t, 0.0), minlength=len(aa))
    denom = np.bincount(row, weights=np.where(use, nu[col] * E, 0.0), minlength=len(aa))
    ok_row = denom > 0
    mu = np.zeros(len(aa))
    mu[ok_row] = num[ok_row] / denom[ok_row]
    row_error = float(np.max(np.abs(rowsum + mu - aa)[ok_row])) if ok_row.any() else 0.0
    entries = ok_row[row] & use
    shape_error = float(np.max(np.abs(t - mu[row] * nu[col] * E)[entries])) if entries.any() else 0.0
    check(row_error < tol and shape_error < tol,
          f'Solver has not converged to the mass constraints (row {row_error:.2e}, form {shape_error:.2e})')
    return dict(row_error=row_error, shape_error=shape_error, rows_checked=int(ok_row.sum()),
                cols_checked=int(ok_col.sum()), rows=len(aa), cols=len(bb))


KKT_LADDER = (1e-7, 1e-9, 1e-11, 1e-13, 1e-15)


def sparse_kkt_check(row, col, c, t, aa, bb, eps, rho, tol=1e-7, shape_tol=1e-8):
    """Optimality check covering every row and column, at any rho/eps (complements sparse_mass_check).

    The solver returns T_ij = exp((f_i + g_j - C_ij)/eps); the unused masses are
    mu_i = exp((f_i - rho)/eps) and nu_j = exp((g_j - rho)/eps). So log T + C/eps = phi_i + psi_j
    with phi = f/eps, psi = g/eps. phi and psi are recovered along a spanning tree of each
    connected component of the positive entries, up to one gauge constant gamma per component
    (phi + gamma, psi - gamma); all other entries test the additive form (shape_error, log scale).
    T is optimal within tol iff some gamma makes every row |T 1 + mu - a| <= tol and every column
    |T^T 1 + nu - b| <= tol. Each condition bounds gamma from one or both sides; the check passes
    if the bounds intersect in every component. Unlike sparse_mass_check, mu and nu need not be
    resolvable in double precision, so no row or column is skipped. Normalized scale, as there.
    Also reports kkt_tol, the smallest tolerance in KKT_LADDER that passes.
    """
    from scipy import sparse
    from scipy.sparse.csgraph import breadth_first_order, connected_components
    row, col, c, t = (np.asarray(v) for v in (row, col, c, t))
    R, K = len(aa), len(bb)
    pos = t > 0
    r, s = row[pos], col[pos]
    L = np.log(t[pos]) + c[pos] / eps
    nodes = R + K
    G = sparse.coo_matrix((np.ones(len(r)), (r, R + s)), shape=(nodes, nodes)).tocsr()
    ncomp, label = connected_components(G, directed=False)
    reps = np.unique(label, return_index=True)[1]
    root = nodes                                     # virtual root joined to one node per component
    G = sparse.coo_matrix((np.ones(len(r) + ncomp), (np.r_[r, np.full(ncomp, root)], np.r_[R + s, reps])),
                          shape=(nodes + 1, nodes + 1)).tocsr()
    order, pred = breadth_first_order(G, root, directed=False, return_predecessors=True)
    order = order[1:].astype(np.int64)               # csgraph returns int32; row*K + col needs int64
    p = pred[order].astype(np.int64)
    tree = p != root
    rnode = np.where(order < R, order, p)[tree]
    cnode = np.where(order < R, p, order)[tree] - R
    keys = r * K + s
    sorter = np.argsort(keys)
    where = sorter[np.searchsorted(keys, rnode * K + cnode, sorter=sorter)]
    check(np.array_equal(keys[where], rnode * K + cnode), 'Spanning-tree edge lookup failed')
    edge_value = np.zeros(len(order))
    edge_value[tree] = L[where]
    value = np.zeros(nodes + 1)
    for v, u, w, is_edge in zip(order.tolist(), p.tolist(), edge_value.tolist(), tree.tolist()):
        value[v] = w - value[u] if is_edge else 0.0
    phi, psi = value[:R], value[R:nodes]
    shape_error = float(np.max(np.abs(L - phi[r] - psi[s]))) if len(r) else 0.0
    check(shape_error < shape_tol, f'Transport is not of the form exp((f + g - C)/eps): log error {shape_error:.3g}')
    delta = aa - np.bincount(row, weights=t, minlength=R)      # = mu at the optimum
    gap = bb - np.bincount(col, weights=t, minlength=K)        # = nu at the optimum
    lab_r, lab_c = label[:R], label[R:]
    m = rho / eps

    def feasible(tol):
        if delta.min() <= -tol or gap.min() <= -tol:
            return False
        lower, upper = np.full(ncomp, -np.inf), np.full(ncomp, np.inf)
        np.minimum.at(upper, lab_r, np.log(delta + tol) - phi + m)            # mu <= delta + tol
        np.maximum.at(lower, lab_c, psi - m - np.log(gap + tol))              # nu <= gap + tol
        hi = delta > tol                                                      # mu >= delta - tol
        np.maximum.at(lower, lab_r[hi], np.log(delta[hi] - tol) - phi[hi] + m)
        hi = gap > tol                                                        # nu >= gap - tol
        np.minimum.at(upper, lab_c[hi], psi[hi] - m - np.log(gap[hi] - tol))
        return bool((lower <= upper).all())

    check(feasible(tol), f'Marginal conditions fail for every gauge at tol {tol:g}: solver has not converged')
    passed = [x for x in KKT_LADDER if feasible(x)]
    return dict(kkt_shape_error=shape_error, kkt_components=int(ncomp), kkt_rows=R, kkt_cols=K,
                kkt_tol=float(min(passed)))


def solve_sparse(row, col, c, a, b, eps=0.2, rho=1.0, nitermax=50000, stopthr=1e-11):
    """Sparse counterpart of solve(): only the stored (row, col) entries are allowed.

    Zero costs must be stored explicitly. Returns transport on the same entries, in the
    original budget units, plus the mass-check diagnostics.
    """
    from scipy import sparse
    a, b = np.asarray(a, float), np.asarray(b, float)
    row, col, c = np.asarray(row, np.int64), np.asarray(col, np.int64), np.asarray(c, float)
    check(a.ndim == b.ndim == 1 and row.shape == col.shape == c.shape, 'Wrong dimensions')
    check(np.isfinite(a).all() and np.isfinite(b).all() and (a >= 0).all() and (b >= 0).all(),
          'Invalid expression budgets')
    check(np.isfinite([eps, rho]).all() and min(eps, rho) > 0, 'Invalid OT parameters')
    check(len(c) == 0 or (np.isfinite(c).all() and (c >= 0).all()), 'Invalid costs')
    check(len(row) == 0 or (row.min() >= 0 and row.max() < len(a) and col.min() >= 0 and col.max() < len(b)),
          'Entry outside the budget vectors')
    check(len(np.unique(row * len(b) + col)) == len(row), 'Duplicate entries')
    T = np.zeros(len(c))
    scale = max(float(a.sum()), float(b.sum()))
    keep = (a[row] > 0) & (b[col] > 0) if len(c) else np.zeros(0, bool)
    if not scale or not keep.any():
        return T, dict(entries=0)
    ia, r = np.unique(row[keep], return_inverse=True)
    ib, s = np.unique(col[keep], return_inverse=True)
    aa, bb = a[ia] / scale, b[ib] / scale
    C = sparse.coo_matrix((c[keep], (r, s)), shape=(len(ia), len(ib)))
    K = backend_sparse()(aa, bb, C, eps, rho, nitermax=nitermax, stopthr=stopthr)
    check(np.array_equal(K.row, C.row) and np.array_equal(K.col, C.col), 'Solver changed the entry order')
    check(np.isfinite(K.data).all(), 'Nonfinite transport')
    info = sparse_mass_check(r, s, c[keep], K.data, aa, bb, eps, rho)
    info.update(sparse_kkt_check(r, s, c[keep], K.data, aa, bb, eps, rho))
    T[keep] = K.data * scale
    info.update(entries=int(keep.sum()), active_rows=len(ia), active_cols=len(ib))
    return T, info


def solve(C, a, b, eps=0.2, rho=1.0):
    a, b, C = np.asarray(a, float), np.asarray(b, float), np.asarray(C, float)
    check(a.ndim == b.ndim == 1 and C.shape == (len(a), len(b)), 'Wrong dimensions')
    check(np.isfinite(a).all() and np.isfinite(b).all() and
          (a >= 0).all() and (b >= 0).all(), 'Invalid expression budgets')
    check(np.isfinite([eps, rho]).all() and min(eps, rho) > 0, 'Invalid OT parameters')
    check(not np.isnan(C).any() and (C >= 0).all(), 'Invalid costs')
    T = np.zeros_like(C)
    scale = max(float(a.sum()), float(b.sum()))
    ia, ib = np.flatnonzero(a > 0), np.flatnonzero(b > 0)
    if not scale or not len(ia) or not len(ib):
        return T
    c = C[np.ix_(ia, ib)]
    ka, kb = np.isfinite(c).any(axis=1), np.isfinite(c).any(axis=0)
    ia, ib = ia[ka], ib[kb]
    if not len(ia) or not len(ib):
        return T
    c = C[np.ix_(ia, ib)]
    aa, bb = a[ia] / scale, b[ib] / scale
    # Equal entropy coefficients for transport AND both unmatched amounts.
    f, g = backend()(aa, bb, c, eps, rho, nitermax=50000,
                     stopthr=1e-11, output_fg=True)
    z = np.exp((f[:, None] + g[None, :] - c) / eps)
    mu, nu = np.exp((f-rho)/eps), np.exp((g-rho)/eps)
    check(np.isfinite(z).all(), 'Nonfinite transport')
    check(np.max(np.abs(z.sum(axis=1) + mu - aa)) < 1e-7 and
          np.max(np.abs(z.sum(axis=0) + nu - bb)) < 1e-7,
          'Solver has not converged to the mass constraints')
    check((z[~np.isfinite(c)] == 0).all(), 'Forbidden transport')
    T[np.ix_(ia, ib)] = z * scale
    return T


def kernel(d, family, reference, r0=None):
    d = np.asarray(d, float)
    check(np.isfinite(reference) and reference > 0 and np.isfinite(d).all()
          and (d >= reference * (1-1e-12)).all(), 'Invalid distances/reference')
    if family == 'inverse':
        return np.minimum(reference / d, 1)
    check(family in ['r0_raw', 'r0_shape'], 'Unknown kernel')
    check(r0 is not None and np.isfinite(r0) and r0 > 0, 'Invalid r0')
    return np.minimum((r0 if family == 'r0_raw' else r0 + reference)/(r0+d), 1)


def contact_neighbors(xy, D, radius, method='radius'):
    """Spatial eligibility; a graph edge does not prove membrane contact."""
    allowed = (D > 0) & (D <= radius)
    if method == 'radius':
        return allowed
    check(method == 'delaunay', 'Unknown contact_graph')
    try:
        triangles = Delaunay(xy).simplices
    except QhullError as exc:
        raise ValueError('Delaunay needs non-collinear 2D coordinates') from exc
    edges = np.zeros(D.shape, dtype=bool)
    for triangle in triangles:
        for i in triangle:
            edges[i, triangle] = True
    return edges & allowed


def cost(D, pairs, ligands, receptors, scores, family, r0, reference, cr, dr,
         contact_edges=None):
    scores = np.asarray(scores, float)
    check(scores.shape == (len(pairs),) and np.isfinite(scores).all() and
          ((scores >= 0) & (scores <= 1)).all(), 'Missing/invalid PPI score')
    n = len(D)
    C = np.full((len(ligands)*n, len(receptors)*n), np.inf)
    for row, s in zip(pairs.itertuples(index=False), scores):
        check(row.type in ['contact', 'diffusion'], 'Unknown LR type')
        eligible = (D > 0) & (D <= (cr if row.type == 'contact' else dr))
        if row.type == 'contact' and contact_edges is not None:
            check(contact_edges.shape == D.shape, 'Contact graph shape mismatch')
            eligible &= contact_edges
        F = np.ones(eligible.sum()) if row.type == 'contact' else kernel(
            D[eligible], family, reference, r0)
        i, j = ligands.index(row.ligand), receptors.index(row.receptor)
        block = C[i*n:(i+1)*n, j*n:(j+1)*n]
        block[eligible] = -np.log(F) - np.log(max(float(s), 1e-8))
    return C


def score_designs(types, s, cfg):
    """Score vectors compared at a fixed distance design.

    no_ppi: S=1. ppi: input scores. ppi_constant (optional): every pair gets the
    mean score of its interaction type, separating score LEVEL from per-pair
    DIFFERENCES. Permutations stay within one type (contact and diffusion costs
    differ structurally). 'random' keeps the v1 behaviour (shuffled_k, both types
    shuffled together); 'exhaustive' enumerates every distinct permutation of one
    type while the other type keeps its input scores (perm_<type>_k).
    """
    types = np.asarray(types)
    s = np.asarray(s, float)
    cases = {'ppi': (s, None)}
    if cfg.get('constant_control', False):
        constant = s.copy()
        for kind in ['contact', 'diffusion']:
            ids = np.flatnonzero(types == kind)
            if len(ids):
                constant[ids] = s[ids].mean()
        cases['ppi_constant'] = (constant, None)
    mode = cfg.get('permutation_mode', 'random')
    check(mode in ['random', 'exhaustive'], 'permutation_mode must be random or exhaustive')
    if mode == 'random':
        rng = np.random.default_rng(cfg.get('seed', 20261112))
        for repeat in range(int(cfg.get('n_permutations', 3))):
            shuffled = s.copy()
            for kind in ['contact', 'diffusion']:
                ids = np.flatnonzero(types == kind)
                shuffled[ids] = rng.permutation(s[ids])
            if not np.array_equal(shuffled, s):
                cases[f'shuffled_{repeat}'] = (shuffled, 'both')
        return cases
    from itertools import permutations
    limit = int(cfg.get('max_permutations', 720))
    for kind in ['contact', 'diffusion']:
        ids = np.flatnonzero(types == kind)
        distinct = sorted(set(permutations(s[ids].tolist())) - {tuple(s[ids].tolist())})
        check(len(distinct) <= limit, f'{len(distinct)} {kind} permutations exceed max_permutations; '
              'use permutation_mode=random')
        for k, values in enumerate(distinct):
            permuted = s.copy()
            permuted[ids] = values
            cases[f'perm_{kind}_{k:03d}'] = (permuted, kind)
    return cases


def numeric_table(path):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        header = next(csv.reader(stream))
    check(len(header) > 1 and all(h.strip() for h in header[1:]) and
          len(set(header[1:])) == len(header[1:]), 'Invalid/duplicate columns')
    x = pd.read_csv(path, index_col=0, dtype=str)
    check(x.index.is_unique and not x.index.hasnans, 'Invalid spot IDs')
    x.index = x.index.astype(str)
    check(all(v.strip() for v in x.index), 'Empty spot ID')
    return x.apply(pd.to_numeric, errors='raise')


def run(folder, output):
    folder, output = Path(folder), Path(output)
    cfg = json.loads((folder/'config.json').read_text())
    check(bool(cfg.get('coordinate_unit')), 'Specify coordinate_unit')
    check(cfg.get('expression_level') in ['gene', 'isoform', 'synthetic'],
          'Specify expression_level')
    x, spatial = numeric_table(folder/'expression.csv'), numeric_table(folder/'spatial.csv')
    check(len(x) >= 2 and set(x.index) == set(spatial.index), 'Spot IDs must match')
    values = x.to_numpy(float)
    check(np.isfinite(values).all() and (values >= 0).all(), 'Invalid expression')
    if cfg.get('expression_state') == 'raw_counts':
        check(cfg.get('full_transcriptome') is True, 'Normalize BEFORE LR subsetting')
        check((values.sum(axis=1) > 0).all(), 'Zero-count spots')
        values = np.log1p(values / values.sum(axis=1)[:, None] * 10000)
        x = pd.DataFrame(values, index=x.index, columns=x.columns)
    else:
        check(cfg.get('expression_state') == 'processed', 'Specify expression_state')
    xy = spatial.loc[x.index, ['x', 'y']].to_numpy(float)
    check(np.isfinite(xy).all() and len(np.unique(xy, axis=0)) == len(xy),
          'Nonfinite/duplicate coordinates')
    lr = pd.read_csv(folder/'lr.csv', dtype=str)[['ligand', 'receptor', 'type']]
    mapping = {'diffusion': 'diffusion', 'contact': 'contact',
               'Secreted Signaling': 'diffusion', 'Cell-Cell Contact': 'contact'}
    custom = cfg.get('type_mapping', {})
    check(not custom or bool(cfg.get('type_mapping_provenance')), 'Document custom type mapping')
    lr['type'] = lr['type'].map({**mapping, **custom})
    check(not lr.isna().any().any() and lr.type.isin(['diffusion', 'contact']).all(),
          'Unknown LR type/ID; ECM annotations require an explicit mapping')
    lr = lr.drop_duplicates()
    check(not lr.duplicated(['ligand', 'receptor']).any(), 'Conflicting LR types')
    positive = set(x.columns[x.sum() > 0])
    keep = lr.ligand.isin(positive) & lr.receptor.isin(positive)
    excluded, lr = lr.loc[~keep], lr.loc[keep].reset_index(drop=True)
    check(len(lr) > 0, 'No expressed LR pairs')
    ligands, receptors = sorted(lr.ligand.unique()), sorted(lr.receptor.unique())
    n = len(x)
    check(len(ligands)*len(receptors)*n*n <= 4000000,
          'Dense prototype too large; use a small documented ROI/LR subset')
    D = cdist(xy, xy)
    nearest = D.copy(); np.fill_diagonal(nearest, np.inf)
    reference, h = float(nearest.min()), float(np.median(nearest.min(axis=1)))
    cr = float(cfg.get('contact_radius_h', 1.1))*h
    diffusion_radius = cfg.get('diffusion_radius_h', 5)
    dr = np.inf if diffusion_radius is None else float(diffusion_radius)*h
    check(np.isfinite(cr) and not np.isnan(dr) and min(cr, dr) > 0, 'Invalid spatial cutoff')
    contact_edges = contact_neighbors(xy, D, cr, cfg.get('contact_graph', 'radius'))
    multiples = sorted(set(cfg.get('r0_h', [0.5, 1, 2, 4])))
    check(len(multiples) > 0 and np.isfinite(multiples).all() and min(multiples) > 0, 'Invalid r0 grid')
    score_cases = {'no_ppi': (np.ones(len(lr)), None)}
    scorefile = folder/'ppi.csv'
    if scorefile.exists():
        check(bool(cfg.get('score_provenance')), 'Specify score_provenance')
        p = pd.read_csv(scorefile, dtype={'ligand': str, 'receptor': str})
        check(not p.duplicated(['ligand', 'receptor']).any(), 'Duplicate PPI rows')
        all_s = p.score.to_numpy(float)
        check(np.isfinite(all_s).all() and ((all_s >= 0) & (all_s <= 1)).all(), 'Invalid scores')
        s = lr.merge(p[['ligand','receptor','score']], on=['ligand','receptor'],
                     how='left', validate='one_to_one').score.to_numpy(float)
        check(np.isfinite(s).all(), 'Missing PPI coverage; no neutral imputation')
        score_cases.update(score_designs(lr.type.to_numpy(), s, cfg))
    backend()  # Verify the pinned backend before creating output files.
    output.mkdir(parents=True, exist_ok=False)
    excluded.to_csv(output/'excluded_lr.csv', index=False)
    lr.to_csv(output/'used_lr.csv', index=False)
    a, b = x[ligands].to_numpy().T.reshape(-1), x[receptors].to_numpy().T.reshape(-1)
    families = cfg.get('r0_families', ['r0_raw', 'r0_shape'])
    check(families and len(set(families)) == len(families)
          and set(families) <= {'r0_raw', 'r0_shape'}, 'Invalid r0 families')
    designs = [(f, m*h) for m in multiples for f in families]
    if cfg.get('include_inverse', True):
        designs = [('inverse', None)] + designs
    baseline, records = None, []
    for family, r0 in designs:
        for score_name, (s, permuted_type) in score_cases.items():
            name = f'case_{len(records):03d}'
            C = cost(D, lr, ligands, receptors, s, family, r0, reference, cr, dr, contact_edges)
            T = solve(C, a, b, cfg.get('eps', 0.2), cfg.get('rho', 1.0))
            if baseline is None:
                baseline = T.copy()
            spot = T.reshape(len(ligands), n, len(receptors), n).sum(axis=(0, 2))
            total = float(T.sum())
            np.savez_compressed(output/f'{name}.npz', C=C, T=T, a=a, b=b,
                                spots=x.index.to_numpy(str), ligands=np.array(ligands),
                                receptors=np.array(receptors), scores=np.asarray(s, float),
                                pair_ligand=lr.ligand.to_numpy(str),
                                pair_receptor=lr.receptor.to_numpy(str),
                                pair_type=lr.type.to_numpy(str))
            pd.DataFrame(spot, index=x.index, columns=x.index).to_csv(output/f'{name}_spots.csv')
            records.append(dict(case=name, family=family, r0=r0, score_mode=score_name,
                permuted_type=permuted_type,
                total_T=total, unused_ligand=float(a.sum()-total),
                unused_receptor=float(b.sum()-total), scores_clipped=int((s < 1e-8).sum()),
                mean_distance=float((spot*D).sum()/total) if total else None,
                relative_L1_change=float(abs(T-baseline).sum()/max(baseline.sum(), 1e-15))))
    pd.DataFrame(records).to_csv(output/'summary.csv', index=False)
    files = ['config.json','expression.csv','spatial.csv','lr.csv'] + (['ppi.csv'] if scorefile.exists() else [])
    if (folder/'preparation.json').is_file():
        files.append('preparation.json')
    manifest = dict(status='complete', config=cfg, h=h, dref=reference,
        solver_url=URL, solver_sha256=SHA256, cccvelo_preprocessing_verified=False,
        contact_graph=cfg.get('contact_graph', 'radius'),
        contact_directed_edges=int(contact_edges.sum()),
        ppi_status='input_scores_used' if scorefile.exists() else 'not_tested_missing_scores',
        score_modes=list(score_cases), permutation_mode=cfg.get('permutation_mode', 'random'),
        input_sha256={f: hashlib.sha256((folder/f).read_bytes()).hexdigest() for f in files})
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2))


def demo(folder):
    folder = Path(folder); folder.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(7); ids = [f's{i}' for i in range(8)]
    pd.DataFrame(rng.poisson(8, (8,4))+1, index=ids,
                 columns=['L1','L2','R1','R2']).to_csv(folder/'expression.csv')
    pd.DataFrame({'x': np.arange(8), 'y': np.zeros(8)}, index=ids).to_csv(folder/'spatial.csv')
    lr = pd.DataFrame(dict(ligand=['L1','L1','L2','L2'], receptor=['R1','R2','R1','R2'],
                          type=['diffusion','diffusion','contact','contact']))
    lr.to_csv(folder/'lr.csv', index=False)
    lr[['ligand','receptor']].assign(score=[0.9,0.2,0.3,0.8]).to_csv(folder/'ppi.csv', index=False)
    (folder/'config.json').write_text(json.dumps(dict(coordinate_unit='synthetic_unit',
        expression_level='synthetic', expression_state='processed',
        score_provenance='SYNTHETIC mock scores, NOT PLM-interact predictions',
        r0_h=[0.5,1,2,4], eps=0.2, rho=1.0, seed=20261112, n_permutations=3), indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--input', default='demo_input')
    p.add_argument('--output', default='results_demo')
    p.add_argument('--demo', action='store_true')
    args = p.parse_args()
    if args.demo:
        demo(args.input)
    run(args.input, args.output)
