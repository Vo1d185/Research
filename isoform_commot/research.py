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


def backend():
    check(TARGET.is_file(), 'Run python fetch_solver.py first')
    check(hashlib.sha256(TARGET.read_bytes()).hexdigest() == SHA256,
          'Official solver checksum mismatch')
    spec = importlib.util.spec_from_file_location('official_unot', TARGET)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.unot_sinkhorn_l1_dense


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
    score_cases = {'no_ppi': np.ones(len(lr))}
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
        score_cases['ppi'] = s
        rng = np.random.default_rng(cfg.get('seed', 20261112))
        for repeat in range(int(cfg.get('n_permutations', 3))):
            shuffled = s.copy()
            for kind in ['contact', 'diffusion']:
                ids = np.flatnonzero(lr.type.to_numpy() == kind)
                shuffled[ids] = rng.permutation(s[ids])
            if not np.array_equal(shuffled, s):
                score_cases[f'shuffled_{repeat}'] = shuffled
    backend()  # Verify the pinned backend before creating output files.
    output.mkdir(parents=True, exist_ok=False)
    excluded.to_csv(output/'excluded_lr.csv', index=False)
    lr.to_csv(output/'used_lr.csv', index=False)
    a, b = x[ligands].to_numpy().T.reshape(-1), x[receptors].to_numpy().T.reshape(-1)
    families = cfg.get('r0_families', ['r0_raw', 'r0_shape'])
    check(families and len(set(families)) == len(families)
          and set(families) <= {'r0_raw', 'r0_shape'}, 'Invalid r0 families')
    designs = [('inverse', None)] + [(f, m*h) for m in multiples for f in families]
    baseline, records = None, []
    for family, r0 in designs:
        for score_name, s in score_cases.items():
            name = f'case_{len(records):03d}'
            C = cost(D, lr, ligands, receptors, s, family, r0, reference, cr, dr, contact_edges)
            T = solve(C, a, b, cfg.get('eps', 0.2), cfg.get('rho', 1.0))
            if baseline is None:
                baseline = T.copy()
            spot = T.reshape(len(ligands), n, len(receptors), n).sum(axis=(0, 2))
            total = float(T.sum())
            np.savez_compressed(output/f'{name}.npz', C=C, T=T, a=a, b=b,
                                spots=x.index.to_numpy(str), ligands=np.array(ligands),
                                receptors=np.array(receptors))
            pd.DataFrame(spot, index=x.index, columns=x.index).to_csv(output/f'{name}_spots.csv')
            records.append(dict(case=name, family=family, r0=r0, score_mode=score_name,
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
