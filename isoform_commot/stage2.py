"""Stage 2: Visium input, sparse collective OT, and dense cross-checks.

Plan: reports/stage2_plan.md (committed before any stage-2 OT run).

  fetch          download a pinned 10x Visium sample (posterior = main, anterior = replicate)
  prepare        COMMOT Visium preprocessing -> input directory (whole slide, or a small region)
  export-genes   full-transcriptome expression of the same spots, for evaluate_targets.py
  run            sparse cost per LR pair + pinned COMMOT unot_sinkhorn_l1_sparse; sparse case files
  compare-dense  rebuild every cost densely from the full distance matrix, solve with the dense
                 solver, and compare cost entries and transport element by element
  pairs          filtered LR pairs of two prepared slides (each slide's own COMMOT filter) and their
                 intersection, with the spot-pair counts at the plan's cutoffs (plan sections 2, 8.3)
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
from research import check, launch_info, solve, solve_sparse

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
COMMOT_PAPER_CUTOFF_UM = 726.1514     # dis_thr = 1000 px in COMMOT's Visium mouse-brain analysis (plan 8.12)
# Main analysis (plan section 8): eps = 0.1, rho = 10, cutoff 363.08 um. Sensitivity: cutoff 500 um;
# eps = 0.2, rho = 1. `run` can override eps, rho and the cutoff; the values used are recorded.
DEFAULTS = dict(eps=0.1, rho=10.0, r0_um=[50.0, 100.0, 200.0, 400.0], diffusion_cutoff_um=COMMOT_TUTORIAL_CUTOFF_UM,
                min_cell_pct=0.05, nitermax=50000, stopthr=1e-11)
OVERRIDES = ['eps', 'rho', 'diffusion_cutoff_um', 'nitermax', 'stopthr']


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


COMMOT_HVG = dict(min_cells=3, min_mean=0.0125, max_mean=3, min_disp=0.5)   # reports/commot_repro_check.md


def commot_hvg(counts, var):
    """COMMOT's downstream-gene notebooks: filter_genes(min_cells=3) on raw counts, log1p, then
    highly_variable_genes(min_mean=0.0125, max_mean=3, min_disp=0.5) (scanpy default flavor 'seurat')."""
    import anndata
    import scanpy as sc
    adata = anndata.AnnData(X=sparse.csr_matrix(counts, dtype=np.float32), var=pd.DataFrame(index=pd.Index(var)))
    sc.pp.filter_genes(adata, min_cells=COMMOT_HVG['min_cells'])
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=COMMOT_HVG['min_mean'], max_mean=COMMOT_HVG['max_mean'],
                                min_disp=COMMOT_HVG['min_disp'])
    return pd.Series(adata.var.highly_variable.reindex(pd.Index(var)).fillna(False).astype(bool).to_numpy(),
                     index=pd.Index(var))


def export_genes(args):
    launch = launch_info()
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
    hvg = commot_hvg(v['X'][idx], v['var'])     # same spots, raw counts
    out.mkdir(parents=True)
    np.savez_compressed(out/'gene_expression.npz', X=X.astype(np.float32), spots=spatial.index.to_numpy(str),
                        genes=np.asarray(v['var'], dtype=str))
    pd.DataFrame(dict(data_symbol=v['var'], lr_symbol=v['var'], detected_fraction=(X > 0).mean(axis=0),
                      commot_hvg=hvg.to_numpy())).to_csv(out/'genes.csv', index=False)
    write_json(out/'gene_export.json', dict(status='complete', launch=launch, input=str(source.resolve()),
        commot_hvg=dict(COMMOT_HVG, n=int(hvg.sum()), flavor='seurat (scanpy default)', on='log1p(raw counts)'),
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
            cfg[key] = int(getattr(args, key)) if key == 'nitermax' else float(getattr(args, key))
    check(cfg['eps'] > 0 and cfg['rho'] > 0 and cfg['diffusion_cutoff_um'] > 0, 'Invalid eps/rho/cutoff')
    return cfg


def run(args):
    launch = launch_info()
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
                   niter=info.get('niter'), nitermax=info.get('nitermax'), cap_reached=info.get('cap_reached'),
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
        print(f'{name}: r0={r0:g} um, {len(T):,} entries, {seconds:.0f}s, {info.get("niter")} iterations', flush=True)
    pd.DataFrame(records).to_csv(out/'summary.csv', index=False)
    from fetch_solver import SHA256, URL
    write_json(out/'manifest.json', dict(status='complete', method='research sparse (stage 2)', launch=launch,
        input=str(source.resolve()), config=cfg,
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


# ------------------------------------------------------------------ pair lists of two slides
def pairs(args):
    """Each slide keeps the pairs that pass its own filter; report both lists and their intersection."""
    from scipy.spatial import cKDTree
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    sources = {'posterior': Path(args.posterior), 'anterior': Path(args.anterior)}
    preps = {k: verified(v) for k, v in sources.items()}
    keys = ['ligand', 'receptor', 'type', 'pathway_name']
    lists = {k: pd.read_csv(v/'lr.csv', dtype=str)[keys] for k, v in sources.items()}
    audits = {k: pd.read_csv(v/'filter_audit.csv', dtype={'ligand': str, 'receptor': str})
              for k, v in sources.items()}
    both = lists['posterior'].merge(lists['anterior'], on=keys, how='outer', indicator=True)
    both['in_posterior'] = both._merge.isin(['both', 'left_only'])
    both['in_anterior'] = both._merge.isin(['both', 'right_only'])
    for k in sources:   # why a pair missing on one slide was excluded there
        reason = audits[k].set_index(['ligand', 'receptor']).reason
        both[f'{k}_filter'] = [reason.get((l, r), 'not_in_database') for l, r in zip(both.ligand, both.receptor)]
    both = both.drop(columns='_merge').sort_values(['type', 'ligand', 'receptor'])
    rows = []
    for kind in ['diffusion', 'contact', 'all']:
        sub_ = both if kind == 'all' else both[both.type == kind]
        p_, a_ = sub_[sub_.in_posterior], sub_[sub_.in_anterior]
        i_ = sub_[sub_.in_posterior & sub_.in_anterior]
        rows.append(dict(type=kind, posterior=len(p_), anterior=len(a_), intersection=len(i_),
                         union=len(sub_), jaccard=len(i_)/len(sub_) if len(sub_) else np.nan,
                         receptors_posterior=p_.receptor.nunique(), receptors_anterior=a_.receptor.nunique(),
                         receptors_intersection=len(set(p_.receptor) & set(a_.receptor))))
    counts = pd.DataFrame(rows)
    edges = []
    for k, v in sources.items():
        sp = pd.read_csv(v/'spatial.csv', index_col=0)
        xy = sp[['x', 'y']].to_numpy(float)
        pr = cKDTree(xy).query_pairs(520.0, output_type='ndarray')
        d = np.linalg.norm(xy[pr[:, 0]] - xy[pr[:, 1]], axis=1)
        for cut in [COMMOT_TUTORIAL_CUTOFF_UM, 363.0, 500.0]:
            edges.append(dict(slide=k, cutoff_um=cut, spot_pairs=int((d <= cut).sum()),
                              max_distance_within=float(d[d <= cut].max()),
                              min_distance_beyond=float(d[d > cut].min()) if (d > cut).any() else np.nan,
                              spots=len(sp), um_per_px=preps[k].get('um_per_px'),
                              neighbour_spacing_um=preps[k].get('neighbour_spacing_um')))
    out.mkdir(parents=True)
    for k in sources:
        both[both[f'in_{k}']][keys].to_csv(out/f'pairs_{k}.csv', index=False)
    both.to_csv(out/'pairs_union.csv', index=False)
    counts.to_csv(out/'counts.csv', index=False)
    pd.DataFrame(edges).to_csv(out/'spot_pairs.csv', index=False)
    write_json(out/'pairs.json', dict(status='complete', launch=launch_info(),
        inputs={k: dict(path=str(v), sample=preps[k]['sample'], prepared_sha256=preps[k]['prepared_sha256'])
                for k, v in sources.items()},
        notes=['Each slide uses only its own filtered pairs (plan: same rules, independent filtering).',
               'spot_pairs counts unordered spot pairs with 0 < d <= cutoff (self pairs excluded).']))
    print(counts.to_string(index=False))
    print(pd.DataFrame(edges).to_string(index=False))


# ------------------------------------------------------------------ robustness (plan section 8.11)
def subsample(args):
    """Keep a random fraction of the spots of a prepared input; values, pairs and coordinates unchanged."""
    launch = launch_info()
    source, out = Path(args.input), Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    check(0 < args.fraction < 1, 'fraction must be in (0, 1)')
    prep = verified(source)
    check('subsample' not in prep, 'Subsample the full prepared input only')
    x = pd.read_csv(source/'expression.csv', index_col=0)
    spatial = pd.read_csv(source/'spatial.csv', index_col=0)
    check(list(x.index) == list(spatial.index), 'Spot order differs between expression and spatial')
    keep = np.sort(np.random.default_rng(args.seed).choice(len(x), size=int(round(args.fraction*len(x))), replace=False))
    out.mkdir(parents=True)
    x.iloc[keep].to_csv(out/'expression.csv')
    spatial.iloc[keep].to_csv(out/'spatial.csv')
    for name in ['lr.csv', 'entities.csv', 'config.json', 'filter_audit.csv']:
        (out/name).write_bytes((source/name).read_bytes())
    files = ['expression.csv', 'spatial.csv', 'lr.csv', 'entities.csv', 'config.json', 'filter_audit.csv']
    meta = {k: v for k, v in prep.items() if k != 'prepared_sha256'}
    meta.update(prepared_sha256={f: digest(out/f) for f in files}, region='random spot subsample',
                spots=int(len(keep)), subsample=dict(parent=str(source.resolve()), parent_prepared_sha256=prep['prepared_sha256'],
                fraction=args.fraction, seed=args.seed, rng='numpy.random.default_rng(seed).choice(n, size, replace=False)',
                launch=launch, note='plan 8.11: values, pairs and micrometre coordinates are not recomputed'))
    write_json(out/'preparation.json', meta)
    print(f'kept {len(keep)} of {len(x)} spots (fraction {args.fraction}, seed {args.seed}) -> {out}')


def robustness(args):
    """Plan 8.11: per receptor, Spearman between subsample and full-data received signal on the kept spots."""
    from scipy.stats import spearmanr
    launch = launch_info()
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    full = {}
    for folder in args.full:
        for case in condition_cases(folder).to_dict('records'):
            full[(case['condition'], case['r0_um'] if pd.notna(case['r0_um']) else None)] = (Path(folder), case['case'])
    rows = []
    for folder in args.subsampled:
        folder = Path(folder)
        manifest = json.loads((folder/'manifest.json').read_text())
        info = json.loads((Path(manifest['input'])/'preparation.json').read_text())['subsample']
        for case in condition_cases(folder).to_dict('records'):
            key = (case['condition'], case['r0_um'] if pd.notna(case['r0_um']) else None)
            check(key in full, f'{folder}: no full-data result for {key}')
            fdir, fcase = full[key]
            sub_plan, full_plan = plans.load(folder/f"{case['case']}.npz"), plans.load(fdir/f'{fcase}.npz')
            sub_pairs = pd.read_csv(folder/'used_lr.csv', dtype=str)
            full_pairs = pd.read_csv(fdir/'used_lr.csv', dtype=str)
            _, sub_rec = plans.received_signal(sub_plan, sub_pairs, len(sub_plan['spots']))
            _, full_rec = plans.received_signal(full_plan, full_pairs, len(full_plan['spots']))
            index = pd.Index(full_plan['spots'].astype(str)).get_indexer(sub_plan['spots'].astype(str))
            check((index >= 0).all(), f'{folder}: subsample spots missing from the full data')
            for receptor in sorted(set(sub_rec) | set(full_rec)):
                s, f = sub_rec.get(receptor), full_rec.get(receptor)
                rho = np.nan
                if s is not None and f is not None and np.ptp(s) > 0 and np.ptp(f[index]) > 0:
                    rho = float(spearmanr(s, f[index]).statistic)
                rows.append(dict(condition=case['condition'], r0_um=case['r0_um'], fraction=info['fraction'],
                                 seed=info['seed'], receptor=receptor, spearman=rho,
                                 missing=('subsample' if s is None else 'full' if f is None else None),
                                 results=str(folder), full_results=str(fdir)))
    table = pd.DataFrame(rows)
    q = lambda p: (lambda s: s.quantile(p))
    summary = table.groupby(['condition', 'r0_um', 'fraction'], dropna=False).agg(
        receptor_runs=('spearman', 'size'), undefined=('spearman', lambda s: int(s.isna().sum())),
        repeats=('seed', 'nunique'), min=('spearman', 'min'), q25=('spearman', q(.25)),
        median=('spearman', 'median'), q75=('spearman', q(.75)), max=('spearman', 'max')).reset_index()
    out.mkdir(parents=True)
    table.to_csv(out/'robustness_receptors.csv', index=False)
    summary.to_csv(out/'robustness_summary.csv', index=False)
    write_json(out/'robustness.json', dict(status='complete', launch=launch, full=args.full, subsampled=args.subsampled,
        rule='plan section 8.11', metric='Spearman(subsample, full) of receptor-level received signal on kept spots'))
    print(summary.to_string(index=False))


CONVERGENCE_COLUMNS = ['condition', 'r0_um', 'niter', 'nitermax', 'cap_reached', 'solves', 'solves_at_cap',
                       'check_passed', 'check_row_error', 'check_shape_error', 'check_kkt_shape_error', 'check_kkt_tol',
                       'check_rows_checked', 'check_rows', 'forbidden_mass', 'max_ligand_budget_excess',
                       'frac_ligand', 'frac_receptor', 'solve_seconds', 'converged', 'results', 'case']


def convergence(args):
    """Plan 8.7 (fourth decision) for result directories, before any target evaluation."""
    launch = launch_info()
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    cases = pd.concat([condition_cases(r) for r in args.results], ignore_index=True)
    check(not cases.duplicated(['condition', 'r0_um'], keep=False).any(), 'A condition/r0 appears twice')
    rows = []
    for case in cases.to_dict('records'):
        case['forbidden_mass'] = forbidden_mass(case)
        case['converged'] = converged(case)
        rows.append({c: case.get(c, np.nan) for c in CONVERGENCE_COLUMNS})
    order = {c: i for i, c in enumerate(CONDITIONS)}
    table = pd.DataFrame(rows).sort_values(['condition', 'r0_um'], key=lambda s: s.map(order) if s.name == 'condition' else s)
    collective = table[table.condition.map(lambda c: CONDITIONS[c][1]) == 'collective']
    out.mkdir(parents=True)
    table.to_csv(out/'convergence.csv', index=False)
    write_json(out/'convergence.json', dict(status='complete', launch=launch, rule='plan section 8.7 (fourth decision)',
        results=args.results, collective_cases=len(collective),
        all_collective_converged=bool(collective.converged.eq(True).all()),
        cases_at_iteration_cap=int(table.cap_reached.eq(True).sum())))
    print(table.drop(columns=['results', 'case']).to_string(index=False))


# ------------------------------------------------------------------ stability (plan section 8.7)
STABILITY = dict(rel_l1=1e-3, spearman=0.999, total=1e-3)


def _cases_by_key(folders):
    out = {}
    for folder in folders:
        for case in condition_cases(folder, strict=False).to_dict('records'):
            key = (case['condition'], None if pd.isna(case['r0_um']) else float(case['r0_um']))
            check(key not in out, f'{folder}: {key} appears twice')
            out[key] = case
    return out


def stability(args):
    """Plan 8.7: 50000- vs 200000-iteration results per receptor; and exact identity of repeated runs."""
    from scipy.stats import spearmanr
    launch = launch_info()
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    short, long_ = _cases_by_key(args.at_50k), _cases_by_key(args.at_200k)
    rows, totals = [], []
    for key, lc in long_.items():
        check(key in short, f'No 50000-iteration result for {key}')
        sc = short[key]
        check(int(sc['nitermax']) == 50000 and int(lc['nitermax']) == 200000, f'{key}: nitermax must be 50000 / 200000')
        loaded = []
        for c in (sc, lc):
            plan = plans.load(Path(c['results'])/f"{c['case']}.npz")
            pairs = pd.read_csv(Path(c['results'])/'used_lr.csv', dtype=str)
            loaded.append((plan, plans.received_signal(plan, pairs, len(plan['spots']))[1]))
        (p50, s50), (p200, s200) = loaded
        check(list(p50['spots']) == list(p200['spots']), f'{key}: spot order differs')
        t50, t200 = float(p50['T'].sum()), float(p200['T'].sum())
        totals.append(dict(condition=key[0], r0_um=key[1], total_T_50k=t50, total_T_200k=t200,
                           rel_total_difference=abs(t50 - t200)/t200, niter_50k=sc.get('niter'), niter_200k=lc.get('niter'),
                           cap_reached_200k=lc.get('cap_reached')))
        for receptor in sorted(set(s50) | set(s200)):
            a, b = s50.get(receptor), s200.get(receptor)
            zero = (a is None or not np.any(a)) and (b is None or not np.any(b))
            l1 = rho = np.nan
            if not zero and a is not None and b is not None:
                l1 = float(np.abs(a - b).sum() / np.abs(b).sum()) if np.abs(b).sum() else np.inf
                if np.ptp(a) > 0 and np.ptp(b) > 0:
                    rho = float(spearmanr(a, b).statistic)
            ok = (not zero) and l1 <= STABILITY['rel_l1'] and rho >= STABILITY['spearman']
            rows.append(dict(condition=key[0], r0_um=key[1], receptor=receptor, both_zero=zero, rel_l1=l1,
                             spearman=rho, meets=None if zero else bool(ok)))
    receptors, totals = pd.DataFrame(rows), pd.DataFrame(totals)
    verdict = []
    for t_ in totals.to_dict('records'):
        r = receptors[(receptors.condition == t_['condition']) & (receptors.r0_um.fillna(-1) == (t_['r0_um'] or -1))]
        used = r[~r.both_zero]
        verdict.append(dict(t_, receptors=len(used), receptors_both_zero=int(r.both_zero.sum()),
            receptors_failing=int((~used.meets.astype(bool)).sum()), max_rel_l1=used.rel_l1.max(),
            min_spearman=used.spearman.min(),
            stable=bool(used.meets.astype(bool).all() and t_['rel_total_difference'] <= STABILITY['total'])))
    verdict = pd.DataFrame(verdict)
    identical = []
    if args.same_a:
        a_, b_ = _cases_by_key(args.same_a), _cases_by_key(args.same_b)
        for key in sorted(set(a_) & set(b_), key=str):
            pa, pb = (plans.load(Path(c['results'])/f"{c['case']}.npz") for c in (a_[key], b_[key]))
            same = all(np.array_equal(pa[k], pb[k]) for k in ['row', 'col', 'T'])
            identical.append(dict(condition=key[0], r0_um=key[1], entries=len(pa['T']), identical=bool(same),
                                  max_abs_T_difference=float(np.abs(pa['T'] - pb['T']).max()) if len(pa['T']) == len(pb['T']) else np.nan))
    identical = pd.DataFrame(identical)
    out.mkdir(parents=True)
    receptors.to_csv(out/'stability_receptors.csv', index=False)
    verdict.to_csv(out/'stability.csv', index=False)
    identical.to_csv(out/'determinism.csv', index=False)
    write_json(out/'stability.json', dict(status='complete', launch=launch, rule='plan section 8.7 (fourth decision)',
        thresholds=STABILITY, at_50k=args.at_50k, at_200k=args.at_200k, same_a=args.same_a, same_b=args.same_b,
        all_stable=bool(verdict.stable.all()), deterministic=bool(identical.identical.all()) if len(identical) else None))
    print(verdict.to_string(index=False)); print(identical.to_string(index=False))


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
    'C726': ('COMMOT', 'mixture', 0.1, 10.0, COMMOT_PAPER_CUTOFF_UM),   # plan 8.12, sensitivity only
}
PAIRED = [('A', 'B'), ('A500', 'B500'), ('Aeps', 'Beps')]   # plan 8.10: same eps, rho, pure collective


def condition_cases(folder, strict=True):
    """One row per case of a result directory, with its plan-section-8 condition code. strict: the solver
    settings must be the plan's (diagnostic runs with another nitermax are read with strict=False)."""
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
    for code in set(codes) if strict else []:              # solver settings as in the plan
        rows = summary[summary.condition == code]
        expected = (10000, 1e-8) if CONDITIONS[code][1] == 'mixture' else (DEFAULTS['nitermax'], DEFAULTS['stopthr'])
        check(np.allclose(rows.nitermax, expected[0]) and np.allclose(rows.stopthr, expected[1]),
              f'{folder}: {code} must use nitermax/stopthr {expected}')
    return summary


def receptor_classes(lr):
    """Plan 8.10: Secreted if a receptor occurs only in diffusion pairs, Contact if only in contact pairs."""
    types = lr.groupby('receptor').type.agg(lambda s: set(s))
    return types.map(lambda s: 'Secreted' if s == {'diffusion'} else 'Contact' if s == {'contact'} else 'both')


def forbidden_mass(case):
    """Plan 8.7: transport on connections the condition forbids (recomputed from the case file)."""
    plan = plans.load(Path(case['results'])/f"{case['case']}.npz")
    if case['method'] == 'ours':
        kinds = plan['pair_type'][plan['pair']]
        allowed = np.where(kinds == 'diffusion', plan['D'] <= case['cutoff_um'] * (1 + 1e-12),
                           plan['D'] <= 1.1 * PITCH_UM)
    else:
        allowed = plan['C'] <= case['cutoff_um'] * (1 + 1e-12)   # C stores the raw distance
    return float(plan['T'][~allowed].sum())


def converged(case):
    """Plan 8.7 (fourth decision): collective conditions pass both mass checks within 1e-7 and carry no
    mass on forbidden connections; reaching the iteration cap is recorded, not a failure. Mixture: n/a."""
    if case['scheme'] != 'collective':
        return np.nan
    kkt = case.get('check_kkt_tol')
    return bool(case.get('check_passed') and pd.notna(kkt) and kkt <= 1e-7 and case.get('forbidden_mass') == 0)


def metric_table(ev, meta, cases, rule):
    k, kind = str(meta['primary_top_tf']), meta['target_score'][0]
    key = 'commot' if rule == 'commot' else k
    sets = pd.read_csv(ev/'target_sets.csv', dtype={'top_tf': str})
    evaluable = set(sets[(sets.top_tf == key) & (sets.status == 'evaluated')].receptor)
    m = pd.read_csv(ev/'case_metrics.csv', dtype={'top_tf': str})
    m = m[(m.level == 'receptor') & (m.top_tf == key) & (m.target_score == kind) & (m.score_mode == 'no_ppi')]
    m = m.merge(cases[['results', 'case', 'condition', 'r0_um']], on=['results', 'case'], how='inner')
    check(set(m.receptor) <= evaluable, 'Metrics for a receptor without an evaluated target set')
    m['r0_key'] = m.r0_um.fillna(0.0)                        # COMMOT cases have no r0
    check(not m.duplicated(['receptor', 'condition', 'r0_key']).any(), 'Duplicate receptor metrics')
    wide = m.set_index(['receptor', 'condition', 'r0_key']).partial_spearman.unstack(['condition', 'r0_key'])
    return m, wide.reindex(sorted(evaluable & set(m.receptor))), evaluable


def select_r0(wide, matrix=False):
    """Plan 8.5: condition A, receptors defined at every r0, highest mean partial Spearman; ties -> smaller r0."""
    a = wide['A']
    a = a.loc[a.notna().all(axis=1)]
    selection = pd.DataFrame(dict(r0_um=a.columns.astype(float), mean_partial_spearman=a.mean().to_numpy(),
                                  receptors=len(a))).sort_values('r0_um').reset_index(drop=True)
    best = selection.mean_partial_spearman.max()
    r0 = float(selection[selection.mean_partial_spearman == best].r0_um.min())
    return (selection, r0, a) if matrix else (selection, r0)


def tolerance_confirmation(a, best, r0, replicates=1000, seed=1):
    """Plan 8.9 (fourth decision): difference of mean partial Spearman between the slide's best r0 and r0,
    and its bootstrap standard error over receptors (both r0 fixed; resampled with replacement)."""
    x, y = a[best].to_numpy(float), a[r0].to_numpy(float)
    idx = np.random.default_rng(seed).integers(0, len(x), size=(replicates, len(x)))
    diffs = x[idx].mean(axis=1) - y[idx].mean(axis=1)
    delta, se = float(x.mean() - y.mean()), float(diffs.std(ddof=1))
    return dict(delta=delta, bootstrap_se=se, replicates=replicates, seed=seed, tolerance_confirmed=bool(delta <= se))


def paired_table(wide, r0, classes):
    """Plan 8.10: per receptor d = ours(r0) - COMMOT; distribution, counts, two-sided exact sign test."""
    from scipy.stats import binomtest
    rows = []
    for x, y in PAIRED:
        if (x, r0) not in wide.columns or y not in wide.columns.get_level_values(0):
            continue
        d_all = (wide[(x, r0)] - wide[y].iloc[:, 0]).dropna()
        for group in ['all', 'Secreted', 'Contact', 'both']:
            d = d_all if group == 'all' else d_all[d_all.index.map(classes) == group]
            pos, neg = int((d > 0).sum()), int((d < 0).sum())
            q = d.quantile([0, .25, .5, .75, 1]).to_numpy() if len(d) else [np.nan]*5
            rows.append(dict(ours=x, baseline=y, r0_um=r0, receptor_class=group, receptors=len(d),
                             ours_higher=pos, baseline_higher=neg, ties=int((d == 0).sum()),
                             sign_test_p=binomtest(pos, pos + neg, 0.5).pvalue if pos + neg else np.nan,
                             mean_difference=d.mean(), min=q[0], q25=q[1], median=q[2], q75=q[3], max=q[4]))
    return pd.DataFrame(rows)


def positive_fraction(m, r0):
    """Plan 8.10: share of receptors with partial Spearman > 0 per condition (A-type: every r0, r0* marked)."""
    rows = []
    for (condition, r0_key), g in m.groupby(['condition', 'r0_key']):
        defined = g.partial_spearman.dropna()
        rows.append(dict(condition=condition, r0_um=r0_key if r0_key else np.nan, selected_r0=bool(r0_key == r0),
                         receptors_defined=len(defined), partial_positive=int((defined > 0).sum()),
                         fraction_partial_positive=float((defined > 0).mean()) if len(defined) else np.nan))
    return pd.DataFrame(rows)


def condition_table(cases, m, wide):
    common = wide.loc[wide.notna().all(axis=1)].index
    rows = []
    for case in cases.to_dict('records'):
        case['forbidden_mass'] = forbidden_mass(case)
        mc = m[(m.results == case['results']) & (m.case == case['case'])]
        mm = mc[mc.receptor.isin(common)]
        g = lambda name: case.get(name, np.nan)
        rows.append(dict(condition=case['condition'], method=case['method'], scheme=case['scheme'], eps=case['eps'],
            rho=case['rho'], cutoff_um=case['cutoff_um'], r0_um=case['r0_um'], frac_ligand=case['frac_ligand'],
            frac_receptor=case['frac_receptor'], total_T=case['total_T'], total_T_contact=g('total_T_contact'),
            total_T_diffusion=g('total_T_diffusion'), mean_distance_um=case['mean_distance'],
            niter=g('niter'), nitermax=g('nitermax'), cap_reached=g('cap_reached'), solves=g('solves'),
            solves_at_cap=g('solves_at_cap'), check_passed=g('check_passed'), kkt_tol=g('check_kkt_tol'),
            rows_checked=g('check_rows_checked'), rows=g('check_rows'),
            max_ligand_budget_excess=g('max_ligand_budget_excess'), forbidden_mass=case['forbidden_mass'],
            row_error=g('check_row_error'), shape_error=g('check_shape_error'), kkt_shape_error=g('check_kkt_shape_error'),
            converged=converged(case),
            solve_seconds=case['solve_seconds'], receptors_defined=int(mc.partial_spearman.notna().sum()),
            receptors_common=len(mm), mean_partial_spearman=mm.partial_spearman.mean(),
            mean_spearman=mm.spearman.mean(), mean_receptor_only_spearman=mm.receptor_only_spearman.mean(),
            results=case['results'], case=case['case']))
    order = {c: i for i, c in enumerate(CONDITIONS)}
    table = pd.DataFrame(rows)
    return table.sort_values(['condition', 'r0_um'], key=lambda s: s.map(order) if s.name == 'condition' else s), common


def report(args):
    """Plan section 8. Posterior: selects r0 (8.5). Anterior (--posterior-report): uses the posterior r0 and
    applies the confirmation rule (8.9). Every table is repeated under COMMOT's target rule (8.8) if evaluated."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from analyze import PALETTE, style_axis
    launch = launch_info()
    ev, out = Path(args.eval), Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    meta = json.loads((ev/'eval_manifest.json').read_text())
    check(meta.get('status') == 'complete', 'Incomplete target evaluation')
    cases = pd.concat([condition_cases(r) for r in meta['results']], ignore_index=True)
    check(not cases.duplicated(['condition', 'r0_um'], keep=False).any(), 'A condition/r0 appears twice')
    lr = pd.read_csv(Path(meta['input'])/'lr.csv', dtype=str)
    classes = receptor_classes(lr)
    m, wide, evaluable = metric_table(ev, meta, cases, 'ours')
    selection, own_r0, a_matrix = select_r0(wide, matrix=True)
    confirmation = None
    if args.posterior_report:                                  # anterior: r0 comes from the posterior (8.9)
        post = json.loads(Path(args.posterior_report).read_text())
        check(post.get('status') == 'complete' and post.get('role') == 'selection', 'Need the posterior report')
        r0_star = float(post['selected_r0_um'])
        ranked = selection.sort_values(['mean_partial_spearman', 'r0_um'], ascending=[False, True]).reset_index(drop=True)
        confirmation = dict(posterior_r0_um=r0_star, rule_selected_r0_here=own_r0, confirmed=bool(own_r0 == r0_star),
            rank_of_posterior_r0=int(ranked.index[ranked.r0_um == r0_star][0]) + 1,
            gap_to_best=float(ranked.mean_partial_spearman.iloc[0] -
                              ranked[ranked.r0_um == r0_star].mean_partial_spearman.iloc[0]),
            receptors=int(selection.receptors.iloc[0]))
        confirmation['strict_confirmed'] = confirmation['confirmed']
        confirmation.update(tolerance_confirmation(a_matrix, own_r0, r0_star))
    else:
        r0_star = own_r0
    table, common = condition_table(cases, m, wide)
    paired = paired_table(wide, r0_star, classes)
    positive = positive_fraction(m, r0_star)
    sets = pd.read_csv(ev/'target_sets.csv', dtype={'top_tf': str})
    if 'rule' not in sets:
        sets['rule'] = 'ours'
    status = sets.groupby(['rule', 'top_tf', 'status']).size().unstack(fill_value=0)
    convergence = table[['condition', 'r0_um', 'niter', 'nitermax', 'cap_reached', 'solves', 'solves_at_cap',
                         'check_passed', 'row_error', 'shape_error', 'kkt_shape_error', 'kkt_tol', 'rows_checked', 'rows',
                         'forbidden_mass', 'max_ligand_budget_excess', 'converged']]
    out.mkdir(parents=True)
    table.to_csv(out/'conditions.csv', index=False)
    convergence.to_csv(out/'convergence.csv', index=False)
    selection.to_csv(out/'r0_selection.csv', index=False)
    paired.to_csv(out/'paired_receptors.csv', index=False)
    positive.to_csv(out/'positive_fraction.csv', index=False)
    status.to_csv(out/'target_set_status.csv')
    classes.rename('receptor_class').to_csv(out/'receptor_classes.csv')
    commot = None
    if (sets.rule == 'commot').any():                           # plan 8.8 sensitivity, same r0
        mc, wc, ec = metric_table(ev, meta, cases, 'commot')
        tc, cc = condition_table(cases, mc, wc)
        tc.to_csv(out/'conditions_commot_rule.csv', index=False)
        paired_table(wc, r0_star, classes).to_csv(out/'paired_receptors_commot_rule.csv', index=False)
        positive_fraction(mc, r0_star).to_csv(out/'positive_fraction_commot_rule.csv', index=False)
        cs = sets[sets.rule == 'commot']
        commot = dict(receptors_evaluable=len(ec), receptors_common=len(cc),
                      complexes_matched=int((cs.complex & (cs.n_linked_tfs > 0)).sum()),
                      complexes_matched_in_commot_order=int((cs.complex & cs.ordered_match).sum()))

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
        if o == 'A' and (sub.r0_um == r0_star).any():
            y = sub[sub.r0_um == r0_star].mean_partial_spearman.iloc[0]
            tag = 'selected' if not args.posterior_report else 'posterior-selected'
            ax.annotate(f'{tag} r0 = {r0_star:g} um', (r0_star, y), xytext=(0, 12), textcoords='offset points',
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
    write_json(out/'report.json', dict(status='complete', launch=launch,
        role='confirmation' if args.posterior_report else 'selection', eval=str(ev), input=meta['input'],
        results=meta['results'], top_tf=meta['primary_top_tf'], target_score=meta['target_score'][0],
        level='receptor', selected_r0_um=r0_star, rule_selected_r0_um=own_r0, r0_rule='plan section 8.5',
        confirmation=confirmation, confirmation_rule='plan section 8.9' if confirmation else None,
        all_converged=bool(table[table.scheme == 'collective'].converged.eq(True).all()),
        convergence_rule='plan section 8.7 (fourth decision)', receptors_evaluable=len(evaluable),
        receptors_common=len(common), receptors_for_r0_selection=int(selection.receptors.iloc[0]),
        receptor_classes=classes.value_counts().to_dict(), commot_rule=commot,
        eval_sha256={f: digest(ev/f) for f in ['case_metrics.csv', 'target_sets.csv']}))
    print(table[['condition', 'r0_um', 'frac_ligand', 'frac_receptor', 'converged', 'mean_partial_spearman']]
          .to_string(index=False))
    print(f'r0 used = {r0_star:g} um (rule here: {own_r0:g}); {len(common)} receptors in every case -> {out}')


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
    r.add_argument('--nitermax', type=int, help='default: config nitermax (diagnostics only; plan: 50000)')
    r.add_argument('--stopthr', type=float, help='default: config stopthr (diagnostics only; plan: 1e-11)')
    c = sub.add_parser('compare-dense'); c.add_argument('--input', required=True); c.add_argument('--results', required=True)
    s = sub.add_parser('pairs'); s.add_argument('--posterior', required=True); s.add_argument('--anterior', required=True)
    s.add_argument('--output', required=True)
    u = sub.add_parser('subsample'); u.add_argument('--input', required=True); u.add_argument('--output', required=True)
    u.add_argument('--fraction', type=float, required=True); u.add_argument('--seed', type=int, required=True)
    b = sub.add_parser('robustness'); b.add_argument('--full', nargs='+', required=True)
    b.add_argument('--subsampled', nargs='+', required=True); b.add_argument('--output', required=True)
    v = sub.add_parser('convergence'); v.add_argument('--results', nargs='+', required=True)
    v.add_argument('--output', required=True)
    y = sub.add_parser('stability'); y.add_argument('--at-50k', nargs='+', required=True)
    y.add_argument('--at-200k', nargs='+', required=True); y.add_argument('--output', required=True)
    y.add_argument('--same-a', nargs='*', default=[], help='determinism check: runs ...')
    y.add_argument('--same-b', nargs='*', default=[], help='... repeated with identical settings')
    o = sub.add_parser('report'); o.add_argument('--eval', required=True, help='evaluate_targets.py output')
    o.add_argument('--output', required=True)
    o.add_argument('--posterior-report', help='anterior only: report.json of the posterior (selected r0; plan 8.9)')
    return p


def main():
    args = parser().parse_args()
    {'fetch': fetch, 'prepare': prepare, 'export-genes': export_genes, 'run': run,
     'compare-dense': compare_dense, 'pairs': pairs, 'subsample': subsample, 'robustness': robustness,
     'stability': stability, 'convergence': convergence,
     'report': report}[args.command](args)


if __name__ == '__main__':
    main()
