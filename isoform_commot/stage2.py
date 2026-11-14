"""Stage 2: Visium input, sparse collective OT, and dense cross-checks.

Plan: reports/stage2_plan.md (committed before any stage-2 OT run).

  fetch          download a pinned 10x Visium sample (posterior = main, anterior = replicate)
  prepare        COMMOT Visium preprocessing -> input directory (whole slide, or a small region)
  export-genes   full-transcriptome expression of the same spots, for evaluate_targets.py
  run            sparse cost per LR pair + pinned COMMOT unot_sinkhorn_l1_sparse; sparse case files
  compare-dense  rebuild every cost densely from the full distance matrix, solve with the dense
                 solver, and compare cost entries and transport element by element
  report         plan section 8: condition table (transported fractions, checks, target metrics),
                 r0 selection (posterior, condition A), paired receptor comparisons, figure

Model (S = 1): diffusion pairs may connect spots with d <= cutoff (d = 0 included),
cost -log(r0/(r0+d)) = log1p(d/r0); contact pairs connect a spot with itself and its
(up to) 6 hexagonal neighbours, cost 0; every other connection is forbidden.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import tarfile
import time

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

import plans
from research import check, solve, solve_sparse

SAMPLES = {
    'posterior': ('V1_Mouse_Brain_Sagittal_Posterior', {
        'filtered_feature_bc_matrix.h5': 'df962209143860488bd679ff7a1b3bb14985208d0929fb814934ba7571c488b6',
        'spatial.tar.gz': 'b878d7eaca487ffc17224a870f86599abf9d45a99033c1ef65661ff119cef0d5'}),
    'anterior': ('V1_Mouse_Brain_Sagittal_Anterior', {
        'filtered_feature_bc_matrix.h5': '56078d8d6fe6c13de248fdb1c518b691cdef78fb00021b659786b4a47c6656d5',
        'spatial.tar.gz': '5f41a803e2bd69fa4dfca6abc8fa2d4e0d76aeb6c72d7038a5fdcf9cc50a36f8'}),
}
URL = 'https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/{sample}/{sample}_{suffix}'
SPOT_DIAMETER_UM = 65.0   # 10x: spot_diameter_fullres spans a theoretical 65 um spot
PITCH_UM = 100.0          # Visium centre-to-centre distance
HEX_OFFSETS = [(0, 2), (0, -2), (1, 1), (1, -1), (-1, 1), (-1, -1)]
CLASSES = {'Secreted Signaling': 'diffusion', 'Cell-Cell Contact': 'contact'}
# COMMOT tutorial dis_thr = 500 px on full-resolution pixels of the posterior slide:
# 500 x 65 um / 89.513 px = 363.0757 um (plan section 8.3); used as a physical distance on both slides.
COMMOT_TUTORIAL_CUTOFF_UM = 363.0757
# Main analysis (plan section 8): eps = 0.1, rho = 10, cutoff 363.08 um. Sensitivity: cutoff 500 um;
# eps = 0.2, rho = 1. `run` can override eps, rho and the cutoff; the values used are recorded.
DEFAULTS = dict(eps=0.1, rho=10.0, r0_um=[50.0, 100.0, 200.0, 400.0], diffusion_cutoff_um=COMMOT_TUTORIAL_CUTOFF_UM,
                min_cell_pct=0.05, nitermax=50000, stopthr=1e-11)
OVERRIDES = ['eps', 'rho', 'diffusion_cutoff_um']


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False, default=str))


# ------------------------------------------------------------------ data
def fetch(args):
    from stage1 import download
    sample, files = SAMPLES[args.sample]
    out = Path(args.output or f'data/visium_mouse_brain_sagittal_{args.sample}')
    for suffix, sha in files.items():
        download(URL.format(sample=sample, suffix=suffix), out/f'{sample}_{suffix}', expected=sha)
    write_json(out/'source.json', dict(status='complete', sample=sample, files=files,
        url=URL.format(sample=sample, suffix='{suffix}'), sha256_source='calculated locally at first download',
        retrieved_utc=datetime.now(timezone.utc).isoformat(timespec='seconds')))
    print(f'{sample} verified -> {out}')


def read_visium(folder, key):
    """Counts (spots x features, CSR), unique var names (anndata make_index_unique), positions, scalefactors."""
    import h5py
    from anndata.utils import make_index_unique
    sample, files = SAMPLES[key]
    folder = Path(folder)
    paths = {s: folder/f'{sample}_{s}' for s in files}
    for suffix, sha in files.items():
        check(paths[suffix].is_file() and digest(paths[suffix]) == sha, f'{paths[suffix]} missing or changed; run fetch')
    with h5py.File(paths['filtered_feature_bc_matrix.h5'], 'r') as f:
        g = f['matrix']
        M = sparse.csc_matrix((g['data'][:], g['indices'][:], g['indptr'][:]), shape=tuple(g['shape'][:]))
        names = [x.decode() for x in g['features']['name'][:]]
        barcodes = [x.decode() for x in g['barcodes'][:]]
    with tarfile.open(paths['spatial.tar.gz']) as tar:
        pos = pd.read_csv(io.BytesIO(tar.extractfile('spatial/tissue_positions_list.csv').read()), header=None,
                          names=['barcode', 'in_tissue', 'array_row', 'array_col', 'pxl_row', 'pxl_col'])
        scalefactors = json.loads(tar.extractfile('spatial/scalefactors_json.json').read())
    pos = pos.set_index('barcode')
    check(pd.Index(barcodes).isin(pos.index).all(), 'Barcodes missing from tissue positions')
    pos = pos.loc[barcodes]
    check((pos.in_tissue == 1).all(), 'Filtered matrix contains spots outside tissue')
    var = pd.Index(make_index_unique(pd.Index(names)))
    return dict(X=M.T.tocsr().astype(np.float64), var=var, barcodes=pd.Index(barcodes), pos=pos,
                scalefactors=scalefactors, sample=sample, n_duplicated=int(pd.Index(names).duplicated().sum()))


def normalize_log1p(X):
    """scanpy normalize_total (default target = median total count) followed by log1p."""
    counts = np.asarray(X.sum(axis=1)).ravel()
    check((counts > 0).all(), 'Spot with zero counts')
    target = float(np.median(counts))
    Y = sparse.diags(target / counts) @ X
    Y = Y.tocsr()
    Y.data = np.log1p(Y.data)
    return Y, target


def hex_pairs(array_row, array_col):
    """Index pairs (i, j), i != j, of spots that are hexagonal neighbours on the Visium array."""
    where = {(r, c): i for i, (r, c) in enumerate(zip(array_row, array_col))}
    check(len(where) == len(array_row), 'Duplicate array positions')
    src, dst = [], []
    for i, (r, c) in enumerate(zip(array_row, array_col)):
        for dr, dc in HEX_OFFSETS:
            j = where.get((r + dr, c + dc))
            if j is not None:
                src.append(i); dst.append(j)
    return np.asarray(src, np.int64), np.asarray(dst, np.int64)


def micrometres(pos, scalefactors):
    um_per_px = SPOT_DIAMETER_UM / float(scalefactors['spot_diameter_fullres'])
    xy = pos[['pxl_col', 'pxl_row']].to_numpy(float) * um_per_px
    s, t = hex_pairs(pos.array_row.to_numpy(), pos.array_col.to_numpy())
    pitch = float(np.median(np.linalg.norm(xy[s] - xy[t], axis=1)))
    check(abs(pitch - PITCH_UM) / PITCH_UM < 0.02,
          f'Neighbour spacing {pitch:.1f} um differs from {PITCH_UM} um by >2%; scale factor not as expected')
    return xy, um_per_px, pitch


def cellchat_pairs(cellchat_dir):
    from cellchat_inventory import load_cellchat
    it, meta = load_cellchat(cellchat_dir)
    cx = pd.read_csv(Path(cellchat_dir)/'complex.csv', dtype=str)
    cols = [c for c in cx.columns if c.startswith('subunit')]
    complexes = {r.complex: [s for s in (getattr(r, c) for c in cols) if isinstance(s, str) and s]
                 for r in cx.itertuples(index=False)}
    pairs = it[it.annotation.isin(CLASSES)].copy()
    pairs['type'] = pairs.annotation.map(CLASSES)
    subunits = {e: complexes.get(e, [e]) for e in set(pairs.ligand) | set(pairs.receptor)}
    return pairs, subunits, meta


def commot_filter(pairs, subunits, var, frac, min_pct):
    """COMMOT filter_lr_database(heteromeric, rule 'min', min_cell_pct): every subunit of both
    entities must be measured and detected in >= min_pct of spots, each subunit separately."""
    pos = {g: i for i, g in enumerate(var)}
    reasons = []
    for row in pairs.itertuples(index=False):
        genes = subunits[row.ligand] + subunits[row.receptor]
        if any(g not in pos for g in genes):
            reasons.append('subunit_not_measured')
        elif any(frac[pos[g]] < min_pct for g in genes):
            reasons.append('subunit_below_min_cell_pct')
        else:
            reasons.append('selected')
    out = pairs.assign(reason=reasons)
    out['min_subunit_fraction'] = [min((frac[pos[g]] for g in subunits[r.ligand] + subunits[r.receptor]
                                        if g in pos), default=0.0) for r in out.itertuples(index=False)]
    return out


def entity_expression(Y, var, entities, subunits):
    pos = {g: i for i, g in enumerate(var)}
    cols = {}
    for e in entities:
        idx = [pos[g] for g in subunits[e]]
        block = Y[:, idx].toarray()
        cols[e] = block.min(axis=1)  # COMMOT heteromeric_rule='min'
    return pd.DataFrame(cols)


def prepare(args):
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    v = read_visium(args.data, args.sample)
    Y, target = normalize_log1p(v['X'])
    frac = np.asarray((v['X'] > 0).sum(axis=0)).ravel() / v['X'].shape[0]   # whole slide
    xy, um_per_px, pitch = micrometres(v['pos'], v['scalefactors'])
    pairs, subunits, cmeta = cellchat_pairs(args.cellchat)
    audit = commot_filter(pairs, subunits, v['var'], frac, args.min_cell_pct)
    selected = audit[audit.reason == 'selected'].copy()
    check(len(selected) and set(selected.type) == {'diffusion', 'contact'}, 'Both LR types must remain')
    spots = np.arange(len(v['barcodes']))
    if args.max_spots:
        centre = np.median(xy, axis=0) if args.center is None else np.asarray(args.center, float)
        spots = np.sort(np.argsort(np.linalg.norm(xy - centre, axis=1), kind='stable')[:args.max_spots])
    if args.pairs_per_type:
        selected = (selected.sort_values(['min_subunit_fraction', 'ligand', 'receptor'], ascending=[False, True, True])
                    .groupby('type', sort=False).head(args.pairs_per_type))
    selected = selected.sort_values(['type', 'ligand', 'receptor']).reset_index(drop=True)
    entities = sorted(set(selected.ligand) | set(selected.receptor))
    expr = entity_expression(Y[spots], v['var'], entities, subunits)
    ids = v['barcodes'][spots]
    expr.index = ids
    pos = v['pos'].iloc[spots]
    spatial = pd.DataFrame(dict(x=xy[spots, 0], y=xy[spots, 1], array_row=pos.array_row.to_numpy(),
                                array_col=pos.array_col.to_numpy()), index=ids)
    nl, nr, n = selected.ligand.nunique(), selected.receptor.nunique(), len(spots)
    out.mkdir(parents=True)
    expr.to_csv(out/'expression.csv', index_label='spot')
    spatial.to_csv(out/'spatial.csv', index_label='spot')
    selected[['ligand', 'receptor', 'type', 'annotation', 'pathway_name', 'interaction_name']].to_csv(
        out/'lr.csv', index=False)
    pd.DataFrame(dict(entity=entities, subunits=[';'.join(subunits[e]) for e in entities])).to_csv(
        out/'entities.csv', index=False)
    audit.drop(columns=['ligand_complex', 'receptor_complex'], errors='ignore')[
        ['ligand', 'receptor', 'type', 'pathway_name', 'reason', 'min_subunit_fraction']].to_csv(
        out/'filter_audit.csv', index=False)
    cfg = dict(DEFAULTS, coordinate_unit='micrometer', expression_level='gene; complexes = subunit minimum',
               expression_state='processed', contact='self + hexagonal neighbours (array_row/array_col)',
               allow_self=True, plan='reports/stage2_plan.md')
    write_json(out/'config.json', cfg)
    files = ['expression.csv', 'spatial.csv', 'lr.csv', 'entities.csv', 'config.json']
    write_json(out/'preparation.json', dict(status='complete', sample=v['sample'], data=str(Path(args.data).resolve()),
        cellchat_commit=cmeta['commit'], cellchat_sha256=cmeta['raw_sha256'],
        normalization=f'normalize_total (target = median total count = {target:g}) then log1p',
        var_names='anndata make_index_unique (first copy keeps the name)', duplicated_features=v['n_duplicated'],
        filter=f'COMMOT filter_lr_database: each subunit detected in >= {args.min_cell_pct:g} of all '
               f'{v["X"].shape[0]} spots', um_per_px=um_per_px,
        spot_diameter_fullres=v['scalefactors']['spot_diameter_fullres'], neighbour_spacing_um=pitch,
        spots=n, region='whole slide' if not args.max_spots else f'{n} spots nearest to {list(np.round(centre, 1))}',
        pairs_selected={k: int(v) for k, v in selected.type.value_counts().items()},
        pairs_per_type_limit=args.pairs_per_type, ligands=int(nl), receptors=int(nr),
        filter_counts={f'{t}/{r}': int(c) for (t, r), c in audit.groupby(['type', 'reason']).size().items()},
        prepared_sha256={f: digest(out/f) for f in files}))
    print(f'Prepared {n} spots, {len(selected)} pairs ({dict(selected.type.value_counts())}), '
          f'{nl} ligands x {nr} receptors -> {out}')


def export_genes(args):
    source, out = Path(args.input), Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    prep = verified(source)
    key = next(k for k, (s, _) in SAMPLES.items() if s == prep['sample'])
    v = read_visium(prep['data'], key)
    Y, _ = normalize_log1p(v['X'])
    spatial = pd.read_csv(source/'spatial.csv', index_col=0)
    idx = v['barcodes'].get_indexer(spatial.index.astype(str))
    check((idx >= 0).all(), 'Spots missing from the Visium matrix')
    X = Y[idx].toarray()
    expression = pd.read_csv(source/'expression.csv', index_col=0)
    _, subunits_of = read_entities(source)
    pos = {g: i for i, g in enumerate(v['var'])}
    for entity in expression.columns:   # re-derive every entity: same spots, normalization and rule
        derived = X[:, [pos[g] for g in subunits_of[entity]]].min(axis=1)
        check(np.allclose(derived, expression[entity].to_numpy(float), rtol=1e-9, atol=1e-12),
              f'{entity} differs from expression.csv')
    out.mkdir(parents=True)
    np.savez_compressed(out/'gene_expression.npz', X=X.astype(np.float32), spots=spatial.index.to_numpy(str),
                        genes=np.asarray(v['var'], dtype=str))
    pd.DataFrame(dict(data_symbol=v['var'], lr_symbol=v['var'], detected_fraction=(X > 0).mean(axis=0))).to_csv(
        out/'genes.csv', index=False)
    write_json(out/'gene_export.json', dict(status='complete', input=str(source.resolve()),
        gene_expression_sha256=digest(out/'gene_expression.npz'), spots=len(idx), genes=len(v['var']),
        note='detected_fraction is computed over these spots; lr_symbol = data symbol (CellChat mouse symbols)'))
    print(f'Exported {len(v["var"])} genes x {len(idx)} spots -> {out}')


def verified(source):
    prep = json.loads((Path(source)/'preparation.json').read_text())
    check(prep.get('status') == 'complete', 'Incomplete preparation')
    for name, sha in prep['prepared_sha256'].items():
        check(digest(Path(source)/name) == sha, f'{name} changed since preparation')
    return prep


def read_entities(source):
    table = pd.read_csv(Path(source)/'entities.csv', dtype=str)
    return table, {r.entity: r.subunits.split(';') for r in table.itertuples(index=False)}


# ------------------------------------------------------------------ sparse model
def spatial_edges(spatial, cutoff):
    """Diffusion edges (d <= cutoff, self included) and contact edges (self + hex neighbours)."""
    xy = spatial[['x', 'y']].to_numpy(float)
    n = len(xy)
    pairs = cKDTree(xy).query_pairs(cutoff, output_type='ndarray')
    s = np.concatenate([np.arange(n), pairs[:, 0], pairs[:, 1]])
    t = np.concatenate([np.arange(n), pairs[:, 1], pairs[:, 0]])
    d = np.linalg.norm(xy[s] - xy[t], axis=1)
    keep = d <= cutoff
    diffusion = (s[keep], t[keep], d[keep])
    hs, ht = hex_pairs(spatial.array_row.to_numpy(), spatial.array_col.to_numpy())
    cs = np.concatenate([np.arange(n), hs]); ct = np.concatenate([np.arange(n), ht])
    contact = (cs, ct, np.linalg.norm(xy[cs] - xy[ct], axis=1))
    return diffusion, contact


def build_sparse_cost(lr, ligands, receptors, a, b, n, edges, r0):
    """Concatenate per-LR-pair sparse blocks; entries with a zero sender or receiver budget are omitted."""
    li, ri = {v: i for i, v in enumerate(ligands)}, {v: i for i, v in enumerate(receptors)}
    rows, cols, costs, dists, kinds = [], [], [], [], []
    for p, pair in enumerate(lr.itertuples(index=False)):
        s, t, d = edges[pair.type]
        r, c = li[pair.ligand]*n + s, ri[pair.receptor]*n + t
        keep = (a[r] > 0) & (b[c] > 0)
        cost = np.log1p(d[keep] / r0) if pair.type == 'diffusion' else np.zeros(int(keep.sum()))
        rows.append(r[keep]); cols.append(c[keep]); costs.append(cost); dists.append(d[keep])
        kinds.append(np.full(int(keep.sum()), p, np.int32))
    cat = lambda x: np.concatenate(x) if x else np.zeros(0)
    return cat(rows).astype(np.int64), cat(cols).astype(np.int64), cat(costs), cat(dists), cat(kinds).astype(np.int64)


def load_input(source):
    source = Path(source)
    cfg = json.loads((source/'config.json').read_text())
    x = pd.read_csv(source/'expression.csv', index_col=0)
    x.index = x.index.astype(str)
    spatial = pd.read_csv(source/'spatial.csv', index_col=0)
    spatial.index = spatial.index.astype(str)
    check(list(x.index) == list(spatial.index), 'Spot order differs between expression and spatial')
    lr = pd.read_csv(source/'lr.csv', dtype=str)
    check(lr.type.isin(['diffusion', 'contact']).all() and not lr.duplicated(['ligand', 'receptor']).any(),
          'Invalid LR table')
    values = x.to_numpy(float)
    check(np.isfinite(values).all() and (values >= 0).all(), 'Invalid expression')
    ligands, receptors = sorted(lr.ligand.unique()), sorted(lr.receptor.unique())
    a, b = x[ligands].to_numpy().T.reshape(-1), x[receptors].to_numpy().T.reshape(-1)
    return cfg, x, spatial, lr, ligands, receptors, a, b


def run_config(cfg, args):
    """Input config with the eps/rho/cutoff given on the command line (if any)."""
    cfg = dict(cfg)
    for key in OVERRIDES:
        if getattr(args, key, None) is not None:
            cfg[key] = float(getattr(args, key))
    check(cfg['eps'] > 0 and cfg['rho'] > 0 and cfg['diffusion_cutoff_um'] > 0, 'Invalid eps/rho/cutoff')
    return cfg


def run(args):
    source, out = Path(args.input), Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    prep = verified(source)
    cfg, x, spatial, lr, ligands, receptors, a, b = load_input(source)
    input_cfg, cfg = cfg, run_config(cfg, args)
    r0s = args.r0_um or cfg['r0_um']
    n = len(x)
    edges = dict(zip(['diffusion', 'contact'], spatial_edges(spatial, cfg['diffusion_cutoff_um'])))
    from research import backend_sparse
    backend_sparse()  # verify the pinned solver before writing anything
    out.mkdir(parents=True)
    lr.to_csv(out/'used_lr.csv', index=False)
    records = []
    for k, r0 in enumerate(r0s):
        name = f'case_{k:03d}'
        row, col, c, d, p = build_sparse_cost(lr, ligands, receptors, a, b, n, edges, r0)
        start = time.time()
        T, info = solve_sparse(row, col, c, a, b, cfg['eps'], cfg['rho'], cfg['nitermax'], cfg['stopthr'])
        seconds = time.time() - start
        types = lr.type.to_numpy()[p]
        plans.save_sparse(out/f'{name}.npz', row, col, c, T, (len(ligands)*n, len(receptors)*n), a, b,
                          x.index, ligands, receptors, D=d, pair=p, pair_ligand=lr.ligand.to_numpy(str),
                          pair_receptor=lr.receptor.to_numpy(str), pair_type=lr.type.to_numpy(str))
        total = float(T.sum())
        rec = dict(case=name, family='r0_raw', r0=float(r0), r0_um=float(r0), score_mode='no_ppi',
                   eps=cfg['eps'], rho=cfg['rho'], diffusion_cutoff_um=cfg['diffusion_cutoff_um'],
                   total_T=total, frac_ligand=total/float(a.sum()), frac_receptor=total/float(b.sum()),
                   unused_ligand=float(a.sum() - total), unused_receptor=float(b.sum() - total),
                   mean_distance=float((T*d).sum()/total) if total else None, entries=len(T),
                   solve_seconds=round(seconds, 1), **{f'check_{k2}': v for k2, v in info.items()})
        for kind in ['diffusion', 'contact']:
            m = types == kind
            rec[f'total_T_{kind}'] = float(T[m].sum())
            rec[f'mean_distance_{kind}'] = float((T[m]*d[m]).sum()/T[m].sum()) if T[m].sum() else None
        records.append(rec)
        print(f'{name}: r0={r0:g} um, {len(T):,} entries, total_T={total:.4g}, {seconds:.0f}s', flush=True)
    pd.DataFrame(records).to_csv(out/'summary.csv', index=False)
    from fetch_solver import SHA256, URL
    write_json(out/'manifest.json', dict(status='complete', method='research sparse (stage 2)', config=cfg,
        overridden={k: input_cfg.get(k) for k in OVERRIDES if input_cfg.get(k) != cfg[k]}, r0_um=r0s,
        spot_pairs=dict(diffusion=int(len(edges['diffusion'][0])), contact=int(len(edges['contact'][0]))),
        h=prep['neighbour_spacing_um'], solver=dict(url=URL, sha256=SHA256,
        function='unot_sinkhorn_l1_sparse'), format='sparse coo (plans.py)', sample=prep['sample'],
        notes=['Entries whose sender or receiver budget is zero are omitted (they cannot carry transport).',
               'Same-spot transport allowed (d = 0).'],
        input_sha256={f: digest(source/f) for f in ['expression.csv', 'spatial.csv', 'lr.csv', 'entities.csv',
                                                    'config.json']}))


# ------------------------------------------------------------------ dense cross-check
def dense_cost(spatial, lr, ligands, receptors, cutoff, r0):
    """Independent construction from the full distance matrix (no KD-tree, no array indices):
    contact = d <= 1.1 x pitch (self + the 6 neighbours on a 100 um hexagonal lattice)."""
    from scipy.spatial.distance import cdist
    xy = spatial[['x', 'y']].to_numpy(float)
    D = cdist(xy, xy)
    n = len(xy)
    C = np.full((len(ligands)*n, len(receptors)*n), np.inf)
    for pair in lr.itertuples(index=False):
        i, j = ligands.index(pair.ligand), receptors.index(pair.receptor)
        block = C[i*n:(i+1)*n, j*n:(j+1)*n]
        if pair.type == 'diffusion':
            ok = D <= cutoff
            block[ok] = np.log1p(D[ok] / r0)
        else:
            block[D <= 1.1 * PITCH_UM] = 0.0
    return C


def compare_dense(args):
    source, results = Path(args.input), Path(args.results)
    verified(source)
    cfg, x, spatial, lr, ligands, receptors, a, b = load_input(source)
    n = len(x)
    check(len(ligands)*len(receptors)*n*n <= 4_000_000, 'Too large for the dense solver; use a small region')
    summary = pd.read_csv(results/'summary.csv')
    cfg = json.loads((results/'manifest.json').read_text())['config']   # values the run actually used
    rows = []
    for case in summary.itertuples(index=False):
        plan = plans.load(results/f'{case.case}.npz')
        C = dense_cost(spatial, lr, ligands, receptors, cfg['diffusion_cutoff_um'], case.r0_um)
        usable = np.isfinite(C) & (a[:, None] > 0) & (b[None, :] > 0)
        dense_entries = set(zip(*np.nonzero(usable)))
        sparse_entries = set(zip(plan['row'].tolist(), plan['col'].tolist()))
        same_pattern = dense_entries == sparse_entries
        cost_diff = float(np.max(np.abs(C[plan['row'], plan['col']] - plan['C']))) if len(plan['C']) else 0.0
        Td = solve(C, a, b, cfg['eps'], cfg['rho'])
        Ts = np.zeros_like(Td)
        Ts[plan['row'], plan['col']] = plan['T']
        rel = float(np.abs(Ts - Td).max() / Td.max()) if Td.max() > 0 else 0.0
        rows.append(dict(case=case.case, r0_um=case.r0_um, entries=len(plan['T']), same_eligible_pattern=same_pattern,
                         max_cost_difference=cost_diff, max_rel_T_difference=rel,
                         total_T_dense=float(Td.sum()), total_T_sparse=float(Ts.sum()),
                         forbidden_mass_dense=float(Td[~np.isfinite(C)].sum()),
                         passed=bool(same_pattern and cost_diff < 1e-12 and rel <= 1e-6)))
    table = pd.DataFrame(rows)
    table.to_csv(results/'dense_comparison.csv', index=False)
    print(table.to_string(index=False))
    check(table.passed.all(), 'Sparse and dense results differ beyond the plan tolerance (1e-6)')


# ------------------------------------------------------------------ report (plan section 8)
CONDITIONS = {   # code: method, scheme, eps, rho, cutoff (um); plan section 8.1
    'A': ('ours', 'collective', 0.1, 10.0, COMMOT_TUTORIAL_CUTOFF_UM),
    'A500': ('ours', 'collective', 0.1, 10.0, 500.0),
    'Aeps': ('ours', 'collective', 0.2, 1.0, COMMOT_TUTORIAL_CUTOFF_UM),
    'B': ('COMMOT', 'collective', 0.1, 10.0, COMMOT_TUTORIAL_CUTOFF_UM),
    'B500': ('COMMOT', 'collective', 0.1, 10.0, 500.0),
    'Beps': ('COMMOT', 'collective', 0.2, 1.0, COMMOT_TUTORIAL_CUTOFF_UM),
    'C': ('COMMOT', 'mixture', 0.1, 10.0, COMMOT_TUTORIAL_CUTOFF_UM),
    'C500': ('COMMOT', 'mixture', 0.1, 10.0, 500.0),
}
PAIRED = [('A', 'B'), ('A', 'C'), ('A500', 'B500'), ('A500', 'C500'), ('Aeps', 'Beps')]


def condition_cases(folder):
    """One row per case of a result directory, with its plan-section-8 condition code."""
    folder = Path(folder)
    manifest = json.loads((folder/'manifest.json').read_text())
    summary = pd.read_csv(folder/'summary.csv')
    ours = manifest.get('method') == 'research sparse (stage 2)'
    if ours:
        summary['cutoff_um'] = summary.diffusion_cutoff_um
        summary['scheme'] = 'collective'
        summary['nitermax'], summary['stopthr'] = manifest['config']['nitermax'], manifest['config']['stopthr']
        summary['check_passed'] = True                     # solve_sparse raises on failure
    else:
        summary['cutoff_um'] = summary.dis_thr_um
        summary['r0_um'] = np.nan
    summary['method'] = 'ours' if ours else 'COMMOT'
    codes = []
    for row in summary.itertuples(index=False):
        hit = [k for k, (m, s, e, r, c) in CONDITIONS.items() if m == row.method and s == row.scheme and
               np.isclose(row.eps, e) and np.isclose(row.rho, r) and np.isclose(row.cutoff_um, c)]
        check(len(hit) == 1, f'{folder}/{row.case}: not a plan-section-8 condition')
        codes.append(hit[0])
    summary['condition'] = codes
    summary['results'] = str(folder)
    for code in set(codes):                                 # solver settings as in the plan
        rows = summary[summary.condition == code]
        expected = (10000, 1e-8) if CONDITIONS[code][1] == 'mixture' else (DEFAULTS['nitermax'], DEFAULTS['stopthr'])
        check(np.allclose(rows.nitermax, expected[0]) and np.allclose(rows.stopthr, expected[1]),
              f'{folder}: {code} must use nitermax/stopthr {expected}')
    return summary


def report(args):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from analyze import PALETTE, style_axis
    ev, out = Path(args.eval), Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    meta = json.loads((ev/'eval_manifest.json').read_text())
    check(meta.get('status') == 'complete', 'Incomplete target evaluation')
    cases = pd.concat([condition_cases(r) for r in meta['results']], ignore_index=True)
    check(not cases.duplicated(['condition', 'r0_um'], keep=False).any(), 'A condition/r0 appears twice')
    k, kind = meta['primary_top_tf'], meta['target_score'][0]
    sets = pd.read_csv(ev/'target_sets.csv')
    evaluable = set(sets[(sets.top_tf == k) & (sets.status == 'evaluated')].receptor)
    m = pd.read_csv(ev/'case_metrics.csv')
    m = m[(m.level == 'receptor') & (m.top_tf == k) & (m.target_score == kind) & (m.score_mode == 'no_ppi')]
    m = m.merge(cases[['results', 'case', 'condition', 'r0_um']], on=['results', 'case'], how='inner')
    check(set(m.receptor) <= evaluable, 'Metrics for a receptor without an evaluated target set')
    m['r0_key'] = m.r0_um.fillna(0.0)                        # COMMOT cases have no r0
    check(not m.duplicated(['receptor', 'condition', 'r0_key']).any(), 'Duplicate receptor metrics')
    wide = m.set_index(['receptor', 'condition', 'r0_key']).partial_spearman.unstack(['condition', 'r0_key'])
    wide = wide.reindex(sorted(evaluable & set(m.receptor)))

    # r0 selection: posterior, condition A, receptors defined at every r0; ties -> smaller r0
    a = wide['A']
    a = a.loc[a.notna().all(axis=1)]
    selection = pd.DataFrame(dict(r0_um=a.columns.astype(float), mean_partial_spearman=a.mean().to_numpy(),
                                  receptors=len(a))).sort_values('r0_um')
    best = selection.mean_partial_spearman.max()
    r0_star = float(selection[selection.mean_partial_spearman == best].r0_um.min())

    # condition table: metrics over the receptors defined in every case
    common = wide.loc[wide.notna().all(axis=1)].index
    rows = []
    for case in cases.itertuples(index=False):
        mc = m[(m.results == case.results) & (m.case == case.case)]
        mm = mc[mc.receptor.isin(common)]
        get = lambda name: getattr(case, name, np.nan)
        rows.append(dict(condition=case.condition, method=case.method, scheme=case.scheme, eps=case.eps,
            rho=case.rho, cutoff_um=case.cutoff_um, r0_um=case.r0_um, frac_ligand=case.frac_ligand,
            frac_receptor=case.frac_receptor, total_T=case.total_T, total_T_contact=get('total_T_contact'),
            total_T_diffusion=get('total_T_diffusion'), mean_distance_um=case.mean_distance,
            check_passed=case.check_passed, kkt_tol=get('check_kkt_tol'), rows_checked=get('check_rows_checked'),
            rows=get('check_rows'), max_ligand_budget_excess=get('max_ligand_budget_excess'),
            solve_seconds=case.solve_seconds, receptors_defined=int(mc.partial_spearman.notna().sum()),
            receptors_common=len(mm), mean_partial_spearman=mm.partial_spearman.mean(),
            mean_spearman=mm.spearman.mean(), mean_receptor_only_spearman=mm.receptor_only_spearman.mean(),
            results=case.results, case=case.case))
    table = pd.DataFrame(rows)
    order = {c: i for i, c in enumerate(CONDITIONS)}
    table = table.sort_values(['condition', 'r0_um'], key=lambda s: s.map(order) if s.name == 'condition' else s)

    # paired receptor comparisons at the selected r0 (descriptive; no p-values)
    paired = []
    for x, y in PAIRED:
        if x not in wide.columns.get_level_values(0) or y not in wide.columns.get_level_values(0):
            continue
        ours_col = wide[(x, r0_star)]
        other = wide[y].iloc[:, 0]
        ok = ours_col.notna() & other.notna()
        d = (ours_col - other)[ok]
        paired.append(dict(ours=x, r0_um=r0_star, baseline=y, receptors=int(ok.sum()),
                           mean_difference=d.mean(), median_difference=d.median(),
                           ours_higher=int((d > 0).sum()), baseline_higher=int((d < 0).sum()), ties=int((d == 0).sum())))
    paired = pd.DataFrame(paired)

    status = sets.groupby(['top_tf', 'status']).size().unstack(fill_value=0)
    out.mkdir(parents=True)
    table.to_csv(out/'conditions.csv', index=False)
    selection.to_csv(out/'r0_selection.csv', index=False)
    paired.to_csv(out/'paired_receptors.csv', index=False)
    status.to_csv(out/'target_set_status.csv')

    # figure: one panel per setting; ours over r0, COMMOT as horizontal references
    panels = [('Main: eps 0.1, rho 10, 363 um', 'A', 'B', 'C'), ('Cutoff 500 um', 'A500', 'B500', 'C500'),
              ('eps 0.2, rho 1, 363 um', 'Aeps', 'Beps', None)]
    panels = [p for p in panels if p[1] in set(table.condition)]
    fig, axes = plt.subplots(1, len(panels), figsize=(3.6*len(panels) + .6, 3.6), sharey=True,
                             constrained_layout=True, squeeze=False)
    for ax, (title, o, b, c) in zip(axes[0], panels):
        sub = table[table.condition == o].sort_values('r0_um')
        ax.plot(sub.r0_um, sub.mean_partial_spearman, color=PALETTE['series1'], lw=2, marker='o', ms=6,
                label='ours (r0)')
        for code, color, style, label in [(b, PALETTE['series2'], '--', 'COMMOT collective'),
                                          (c, PALETTE['series3'], ':', 'COMMOT mixture')]:
            if code and code in set(table.condition):
                ax.axhline(table[table.condition == code].mean_partial_spearman.iloc[0], color=color, ls=style,
                           lw=2, label=label)
        if o == 'A':
            y = sub[sub.r0_um == r0_star].mean_partial_spearman.iloc[0]
            ax.annotate(f'selected r0 = {r0_star:g} um', (r0_star, y), xytext=(0, 12), textcoords='offset points',
                        ha='center', fontsize=8, color=PALETTE['text_secondary'])
        ax.set_xscale('log', base=2)
        ax.set_xticks(sorted(sub.r0_um), [f'{v:g}' for v in sorted(sub.r0_um)])
        ax.minorticks_off()
        ax.set_xlabel('r0 (um)', fontsize=9, color=PALETTE['text_secondary'])
        ax.set_title(title, fontsize=9, loc='left')
        style_axis(ax)
    axes[0, 0].set_ylabel(f'Mean partial Spearman ({len(common)} receptors)', fontsize=9,
                          color=PALETTE['text_secondary'])
    axes[0, 0].legend(fontsize=8, frameon=False, loc='best')
    fig.savefig(out/'partial_spearman.png', dpi=180)
    plt.close(fig)
    write_json(out/'report.json', dict(status='complete', eval=str(ev), results=meta['results'], top_tf=int(k),
        target_score=kind, level='receptor', selected_r0_um=r0_star, r0_rule='plan section 8.5',
        receptors_evaluable=len(evaluable), receptors_common=len(common),
        receptors_for_r0_selection=int(len(a)), eval_sha256={f: digest(ev/f) for f in
                                                            ['case_metrics.csv', 'target_sets.csv']}))
    print(table[['condition', 'r0_um', 'frac_ligand', 'frac_receptor', 'mean_partial_spearman']].to_string(index=False))
    print(f'selected r0 = {r0_star:g} um; {len(common)} receptors in every case -> {out}')


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fetch'); f.add_argument('--sample', choices=list(SAMPLES), required=True)
    f.add_argument('--output')
    q = sub.add_parser('prepare')
    q.add_argument('--sample', choices=list(SAMPLES), default='posterior')
    q.add_argument('--data', required=True); q.add_argument('--cellchat', default='data/cellchatdb')
    q.add_argument('--output', required=True)
    q.add_argument('--min-cell-pct', type=float, default=DEFAULTS['min_cell_pct'])
    q.add_argument('--max-spots', type=int, help='small region: this many spots nearest to the centre')
    q.add_argument('--center', type=float, nargs=2, help='region centre in um (default: median spot position)')
    q.add_argument('--pairs-per-type', type=int, help='small region only: keep this many pairs per type')
    e = sub.add_parser('export-genes'); e.add_argument('--input', required=True); e.add_argument('--output', required=True)
    r = sub.add_parser('run'); r.add_argument('--input', required=True); r.add_argument('--output', required=True)
    r.add_argument('--r0-um', type=float, nargs='+', help='default: config r0_um')
    r.add_argument('--eps', type=float, help='default: config eps')
    r.add_argument('--rho', type=float, help='default: config rho')
    r.add_argument('--diffusion-cutoff-um', type=float, help='default: config diffusion_cutoff_um')
    c = sub.add_parser('compare-dense'); c.add_argument('--input', required=True); c.add_argument('--results', required=True)
    o = sub.add_parser('report'); o.add_argument('--eval', required=True, help='evaluate_targets.py output')
    o.add_argument('--output', required=True)
    return p


def main():
    args = parser().parse_args()
    {'fetch': fetch, 'prepare': prepare, 'export-genes': export_genes, 'run': run,
     'compare-dense': compare_dense, 'report': report}[args.command](args)


if __name__ == '__main__':
    main()
