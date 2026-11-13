"""Pinned CCCvelo LR priors + real h5ad input for a small parameter pilot.

This adapts inputs for the proposed cost; it does not reproduce CCCvelo's
SCT/ALRA, TF/TG network, dynamics, or isoform measurements.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist

from research import check, contact_neighbors, run


CCC_COMMIT = 'fac213c600e7929b5d14e7d95dbdaabd722f9bc9'
LR_FILES = {
    'diff_LigRecDB.csv': ('diffusion', '329db978c93a8edda011acf96fc6c90d192c9012b664b3dbc1d652ac681e79ef'),
    'cont_LigRecDB.csv': ('contact', 'dfc350e2798e46b23a2503d58ea0c0095c91b4a96ca99685f125df665617c1fc'),
}
# Verified against the original provider's machine-readable source, not OCR:
# https://github.com/aristoteleo/spateo-release/blob/1bd8a35e6e7b0dfc712d76b1c32987dc06652774/spateo/sample_data.py#L157
# No author-published checksum was found for this h5ad.
BRAIN_URL = 'https://www.dropbox.com/s/c5tu4drxda01m0u/mousebrain_bin60.h5ad?dl=1'
BRAIN_SOURCE = 'https://github.com/aristoteleo/spateo-release/blob/1bd8a35e6e7b0dfc712d76b1c32987dc06652774/spateo/sample_data.py#L157'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False))


def download(url, target, expected=None, hdf5=False):
    target = Path(target)
    if target.exists():
        check(expected is not None and digest(target) == expected,
              f'Refusing to overwrite {target}; use another path or inspect the existing file')
        return
    target.parent.mkdir(parents=True, exist_ok=True)
    part = target.with_name(target.name + '.part')
    check(not part.exists(), f'Incomplete download exists: {part}; inspect/remove it before retrying')
    try:
        with urlopen(url, timeout=60) as response, part.open('xb') as stream:
            first = response.read(8)
            check(not hdf5 or first == b'\x89HDF\r\n\x1a\n',
                  f'The URL did not return an HDF5 file '
                  f'(HTTP {response.getcode()}, Content-Type: {response.headers.get("Content-Type", "unknown")}). '
                  'No data file was saved. Check the exact source URL and network access.')
            stream.write(first)
            received, next_progress = len(first), 64 * 1024 * 1024
            for block in iter(lambda: response.read(1024 * 1024), b''):
                stream.write(block)
                received += len(block)
                if hdf5 and received >= next_progress:
                    print(f'Downloaded {received / 1024**2:.0f} MiB...', flush=True)
                    next_progress += 64 * 1024 * 1024
            expected_size = response.headers.get('Content-Length')
            check(expected_size is None or received == int(expected_size),
                  'Incomplete download: response size does not match Content-Length')
        check(expected is None or digest(part) == expected, 'Download checksum mismatch')
        if hdf5:
            import h5py
            try:
                with h5py.File(part, 'r') as handle:
                    check({'obs', 'var'} <= set(handle), 'Downloaded HDF5 is not an h5ad expression file')
            except OSError as exc:
                raise ValueError('Downloaded HDF5 is incomplete or unreadable; no final file was saved') from exc
        part.rename(target)
    except Exception:
        part.unlink(missing_ok=True)
        raise


def fetch_lr(folder):
    folder = Path(folder)
    for filename, (_, sha) in LR_FILES.items():
        url = f'https://raw.githubusercontent.com/SunXQlab/CCCvelo/{CCC_COMMIT}/Databases/{filename}'
        download(url, folder / filename, expected=sha)
    print(f'CCCvelo LR files verified at {CCC_COMMIT}: {folder}')


def load_lr(folder):
    frames = []
    for filename, (kind, sha) in LR_FILES.items():
        path = Path(folder) / filename
        check(path.is_file() and digest(path) == sha, f'Wrong/missing pinned LR file: {path}')
        frame = pd.read_csv(path, dtype=str).rename(columns={'source': 'ligand', 'target': 'receptor'})
        frame['type'] = kind
        frames.append(frame[['ligand', 'receptor', 'type', 'origin']])
    pairs = pd.concat(frames, ignore_index=True)
    check(not pairs.isna().any().any(), 'Missing LR identifiers')
    check(not pairs.duplicated(['ligand', 'receptor']).any(), 'Duplicate/conflicting LR pairs')
    return pairs


def open_data(path):
    import anndata as ad
    check(Path(path).is_file(), f'Data file not found: {Path(path).resolve()}. '
          'First run python stage1.py download-data, or pass --h5ad with the existing file path.')
    return ad.read_h5ad(path, backed='r')


def matrix_and_genes(data, source, gene_key=None):
    if source == 'raw':
        check(data.raw is not None, 'No adata.raw matrix')
        matrix, var = data.raw.X, data.raw.var
    elif source == 'X':
        matrix, var = data.X, data.var
    else:
        check(source.startswith('layer:') and source[6:] in data.layers, 'Unknown matrix source')
        matrix, var = data.layers[source[6:]], data.var
    check(gene_key is None or gene_key in var, 'Unknown gene-key')
    genes = var.index.astype(str) if gene_key is None else var[gene_key].astype(str)
    check(not pd.Index(genes).duplicated().any(), 'Duplicate gene identifiers; resolve upstream')
    return matrix, np.asarray(genes, dtype=str)


def dense_rows(matrix, indices):
    x = matrix[np.asarray(indices), :]
    return x.toarray().astype(float) if sparse.issparse(x) else np.asarray(x, dtype=float)


def inspect_data(path, output=None):
    data = open_data(path)
    try:
        matrices = {}
        for source in ['X'] + (['raw'] if data.raw is not None else []) + ['layer:' + k for k in data.layers]:
            m, genes = matrix_and_genes(data, source)
            sample = dense_rows(m, np.arange(min(32, data.n_obs)))
            matrices[source] = dict(shape=list(m.shape), gene_examples=genes[:8].tolist(),
                sampled_min=float(sample.min()), sampled_max=float(sample.max()),
                sampled_integer_like=bool(np.allclose(sample, np.rint(sample), atol=1e-6, rtol=0)),
                sampled_nonnegative=bool(np.isfinite(sample).all() and (sample >= 0).all()))
        report = dict(path=str(Path(path).resolve()), n_obs=data.n_obs, n_vars=data.n_vars,
            matrices=matrices, obs_columns=list(data.obs.columns), var_columns=list(data.var.columns),
            spatial_keys={k: list(data.obsm[k].shape) for k in data.obsm if 'spatial' in k.lower()},
            annotation_examples={k: data.obs[k].astype(str).value_counts().head(15).to_dict()
                for k in ['scc_anno', 'Celltype', 'cell_type', 'Cluster'] if k in data.obs},
            warning='Integer-like samples do NOT establish raw-count provenance. Check source documentation; never normalize processed data twice.')
        if output:
            check(not Path(output).exists(), 'Inspection output exists')
            write_json(output, report)
        print(json.dumps(report, indent=2, ensure_ascii=False))
        return report
    finally:
        data.file.close()


def map_genes(genes, policy, mapping_file=None):
    if mapping_file:
        table = pd.read_csv(mapping_file, dtype=str)
        check({'input_gene', 'lr_gene'} <= set(table), 'Mapping needs input_gene,lr_gene')
        check(not table[['input_gene', 'lr_gene']].isna().any().any()
              and not table.input_gene.duplicated().any() and not table.lr_gene.duplicated().any(),
              'Use an explicit one-to-one gene/ortholog map')
        mapping = dict(zip(table.input_gene, table.lr_gene))
        mapped = [mapping.get(g) for g in genes]
    else:
        mapped = [g.upper() if policy == 'author-uppercase' else g for g in genes]
    valid = [g for g in mapped if g is not None]
    check(len(valid) == len(set(valid)), 'Gene mapping creates duplicate entities; resolve upstream')
    return mapped


def prepare(args):
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    check(args.max_spots >= 4 and args.pairs_per_type >= 1, 'Invalid pilot size')
    check(0 <= args.min_fraction <= 1 and args.min_mean >= 0, 'Invalid expression filter')
    check(np.isfinite(args.contact_radius_h) and args.contact_radius_h > 0, 'Invalid contact radius')
    pairs = load_lr(args.lr_dir)
    data = open_data(args.h5ad)
    try:
        check(data.obs_names.is_unique, 'Duplicate spot IDs')
        check(args.spatial_key in data.obsm, 'Unknown spatial-key; run inspect first')
        xy_all = np.asarray(data.obsm[args.spatial_key], dtype=float)
        check(xy_all.shape == (data.n_obs, 2) and np.isfinite(xy_all).all(), 'Need finite 2D coordinates')
        if args.center:
            center = np.asarray(args.center)
        else:
            center = xy_all[np.argmin(np.linalg.norm(xy_all - np.median(xy_all, axis=0), axis=1))]
        # Compact neighborhood, not a random subsample; keep all annotated types.
        ids = np.sort(np.argsort(np.linalg.norm(xy_all-center, axis=1), kind='stable')[:args.max_spots])
        xy = xy_all[ids]
        check(len(np.unique(xy, axis=0)) == len(xy), 'Duplicate coordinates in ROI')
        matrix, genes = matrix_and_genes(data, args.matrix, args.gene_key)
        values = dense_rows(matrix, ids)
        check(np.isfinite(values).all() and (values >= 0).all(), 'Expression must be finite and nonnegative')
        libraries = values.sum(axis=1)
        check((libraries > 0).all(), 'ROI contains zero-expression spots')
        if args.expression_state == 'counts':
            check(np.allclose(values, np.rint(values), atol=1e-6, rtol=0),
                  'Non-integer values declared counts; inspect expression layers before proceeding')
            # Full library size, including genes absent from the mapping/LR panel.
            values = np.log1p(values / libraries[:, None] * 10000)
        mapped = map_genes(genes, args.symbol_policy, args.gene_map)
        positions = {g: i for i, g in enumerate(mapped) if g is not None}
        fraction, means = (values > 0).mean(axis=0), values.mean(axis=0)
        audit = []
        for row in pairs.itertuples(index=False):
            reason = 'eligible'
            if any(any(c in g for c in [' ', ':', '+', '/', '_']) for g in [row.ligand, row.receptor]):
                reason = 'unsupported_complex_or_nonstandard_entity'
            elif row.ligand not in positions or row.receptor not in positions:
                reason = 'gene_not_mapped_or_measured'
            elif any(fraction[positions[g]] < args.min_fraction or means[positions[g]] < args.min_mean
                     or means[positions[g]] <= 0 for g in [row.ligand, row.receptor]):
                reason = 'below_expression_filter'
            support = min(fraction[positions[g]] for g in [row.ligand, row.receptor]) if reason == 'eligible' else 0.
            audit.append(dict(ligand=row.ligand, receptor=row.receptor, type=row.type,
                              origin=row.origin, reason=reason, expression_support=float(support)))
        audit = pd.DataFrame(audit)
        candidates = audit[audit.reason == 'eligible'].sort_values(
            ['expression_support', 'ligand', 'receptor'], ascending=[False, True, True])
        selected = candidates.groupby('type', sort=False).head(args.pairs_per_type).copy()
        if set(selected.type) != {'diffusion', 'contact'}:
            out.mkdir(parents=True)
            audit.to_csv(out/'lr_audit.csv', index=False)
            write_json(out/'preparation.json', dict(status='blocked',
                reason='Both LR types must be expressed/mapped in the pilot ROI',
                candidates=candidates.groupby('type').size().to_dict()))
            raise ValueError(f'Insufficient LR type coverage. Inspect {out}/lr_audit.csv; '
                             'change ROI/mapping and use a NEW output directory. Do not invent labels.')
        used_genes = sorted(set(selected.ligand) | set(selected.receptor))
        nl, nr, n = selected.ligand.nunique(), selected.receptor.nunique(), len(ids)
        check(nl*nr*n*n <= 4000000, 'Pilot exceeds dense solver capacity; reduce max-spots/pairs-per-type')
        D = cdist(xy, xy); nearest = D.copy(); np.fill_diagonal(nearest, np.inf)
        h = float(np.median(nearest.min(axis=1)))
        edges = contact_neighbors(xy, D, args.contact_radius_h*h, 'delaunay')
        check(edges.any(), 'No contact edges at this cutoff')
        spot_ids = data.obs_names[ids].astype(str)
        expression = pd.DataFrame(values[:, [positions[g] for g in used_genes]],
                                  index=spot_ids, columns=used_genes)
        # Validate optional real scores before writing an input directory.
        ppi = None
        if args.ppi:
            check(bool(args.score_provenance), 'Real PPI scores require --score-provenance')
            scores = pd.read_csv(args.ppi, dtype={'ligand': str, 'receptor': str})
            check(not scores.duplicated(['ligand', 'receptor']).any(), 'Duplicate PPI pairs')
            ppi = selected[['ligand', 'receptor']].merge(scores[['ligand', 'receptor', 'score']],
                on=['ligand', 'receptor'], how='left', validate='one_to_one')
            check(np.isfinite(ppi.score).all() and ppi.score.between(0, 1).all(), 'Missing/invalid PPI coverage')
        cfg = dict(coordinate_unit=args.coordinate_unit, expression_level='gene',
            expression_state='processed', contact_graph='delaunay',
            contact_radius_h=args.contact_radius_h, diffusion_radius_h=None,
            r0_h=[0.5, 1, 2, 4], r0_families=['r0_raw'], eps=0.2, rho=1.0,
            seed=20261113, n_permutations=3,
            lr_source=f'CCCvelo/{CCC_COMMIT}', score_provenance=args.score_provenance)
        selected_keys = set(zip(selected.ligand, selected.receptor))
        audit.loc[(audit.reason == 'eligible') & ~pd.Series(
            [(l, r) in selected_keys for l, r in zip(audit.ligand, audit.receptor)]), 'reason'] = 'outside_fixed_pilot_panel'
        audit.loc[[(l, r) in selected_keys for l, r in zip(audit.ligand, audit.receptor)], 'reason'] = 'selected'
        report = dict(status='complete', stage='CCCvelo-classified LR parameter pilot',
            h5ad=str(Path(args.h5ad).resolve()), h5ad_sha256=digest(args.h5ad),
            cccvelo_commit=CCC_COMMIT, lr_sha256={f: sha for f, (_, sha) in LR_FILES.items()},
            matrix=args.matrix, expression_state=args.expression_state,
            expression_transform='full-library log1p(10000*x/total)' if args.expression_state == 'counts' else 'none',
            gene_key=args.gene_key, symbol_policy=args.symbol_policy,
            gene_map_sha256=digest(args.gene_map) if args.gene_map else None,
            spatial_key=args.spatial_key, coordinate_unit=args.coordinate_unit,
            source_spots=data.n_obs, selected_spots=n, roi_center=center.tolist(),
            roi_radius=float(np.linalg.norm(xy-center, axis=1).max()), h=h,
            filter_min_fraction=args.min_fraction, filter_min_mean=args.min_mean,
            lr_candidates=candidates.groupby('type').size().to_dict(),
            lr_selected=selected.groupby('type').size().to_dict(),
            shared_ligands=int((selected.groupby('ligand').size() > 1).sum()),
            shared_receptors=int((selected.groupby('receptor').size() > 1).sum()),
            contact_directed_edges=int(edges.sum()), contact_isolated_spots=int((~edges.any(axis=1)).sum()),
            ppi_status='provided_scores' if ppi is not None else 'not_tested_missing_scores',
            limitations=['Not a full CCCvelo preprocessing reproduction; no SCT/ALRA or TF/TG filtering.',
                'Author-uppercase is symbol matching, NOT a validated ortholog mapping.',
                'ROI truncates outside senders; Delaunay contacts are spatial proxies.',
                'contact_radius_h is a pilot choice, not the paper physical-distance threshold.',
                'Gene-level input does not measure isoform abundance or interaction frequency.'])
        out.mkdir(parents=True)
        expression.to_csv(out/'expression.csv', index_label='spot')
        pd.DataFrame(xy, index=spot_ids, columns=['x', 'y']).to_csv(out/'spatial.csv', index_label='spot')
        data.obs.iloc[ids].to_csv(out/'cell_meta.csv', index_label='spot')
        selected[['ligand', 'receptor', 'type']].to_csv(out/'lr.csv', index=False)
        audit.to_csv(out/'lr_audit.csv', index=False)
        if ppi is not None:
            ppi.to_csv(out/'ppi.csv', index=False)
        write_json(out/'config.json', cfg)
        report['prepared_sha256'] = {f: digest(out/f) for f in
            ['expression.csv', 'spatial.csv', 'lr.csv'] + (['ppi.csv'] if ppi is not None else [])}
        write_json(out/'preparation.json', report)
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7, 6))
        ax.scatter(xy_all[:, 0], xy_all[:, 1], s=1, color='lightgray', rasterized=True)
        ax.scatter(xy[:, 0], xy[:, 1], s=8, color='#2166ac')
        ax.set_aspect('equal'); ax.set_title('Pilot ROI (blue); outside senders excluded')
        ax.set_xlabel(args.coordinate_unit); ax.set_ylabel(args.coordinate_unit)
        fig.savefig(out/'roi.png', dpi=160); plt.close(fig)
        print(f'Prepared {n} spots, {len(selected)} LR pairs; inspect {out}/roi.png and lr_audit.csv')
    finally:
        data.file.close()


def run_stage(args):
    from analyze import analyze, paired_metrics
    source, out = Path(args.input), Path(args.output)
    check(not out.exists(), 'Output exists')
    preparation = json.loads((source/'preparation.json').read_text())
    check(preparation['status'] == 'complete' and preparation['cccvelo_commit'] == CCC_COMMIT,
          'Expected a prepared stage-one input')
    for name, sha in preparation['prepared_sha256'].items():
        check(Path(name).name == name and digest(source/name) == sha,
              f'Prepared input changed: {name}; prepare a new experiment input')
    check((source/'ppi.csv').exists() == ('ppi.csv' in preparation['prepared_sha256']),
          'PPI input changed since preparation; prepare a new input with score provenance')
    cfg = json.loads((source/'config.json').read_text())
    settings = [(0.2, 1.0)]
    if args.sensitivity:
        settings += [(e, p) for e in [0.1, 0.2, 0.5] for p in [0.5, 1.0, 2.0] if (e, p) not in settings]
    out.mkdir(parents=True)
    all_rows, all_types, stability = [], [], []
    references = {}
    for index, (eps, rho) in enumerate(settings):
        folder = out/f'ot_{index:02d}'
        inputs = folder/'input'
        shutil.copytree(source, inputs)
        current = dict(cfg, eps=eps, rho=rho)
        write_json(inputs/'config.json', current)
        print(f'Running eps={eps}, rho={rho}', flush=True)
        run(inputs, folder/'results')
        analyze(folder/'results', inputs, folder/'analysis')
        summary = pd.read_csv(folder/'results'/'summary.csv')
        for row in summary.itertuples(index=False):
            key = (row.family, None if pd.isna(row.r0) else row.r0, row.score_mode)
            path = folder/'results'/f'{row.case}.npz'
            if index == 0:
                references[key] = path
            with np.load(references[key]) as base, np.load(path) as trial:
                metrics = paired_metrics(base['T'], trial['T'])
            stability.append(dict(eps=eps, rho=rho, family=row.family, r0=row.r0,
                score_mode=row.score_mode, relative_L1_vs_eps02_rho1=metrics['paired_relative_L1'],
                allocation_TV_vs_eps02_rho1=metrics['allocation_TV']))
        all_rows.append(summary.assign(eps=eps, rho=rho, experiment=folder.name))
        all_types.append(pd.read_csv(folder/'analysis'/'by_type.csv').assign(eps=eps, rho=rho))
    pd.concat(all_rows).to_csv(out/'all_conditions.csv', index=False)
    pd.concat(all_types).to_csv(out/'all_types.csv', index=False)
    pd.DataFrame(stability).to_csv(out/'ot_parameter_stability.csv', index=False)
    write_json(out/'stage1_manifest.json', dict(status='complete', experiments=len(settings),
        cccvelo_commit=CCC_COMMIT, input_preparation_sha256=digest(source/'preparation.json'),
        ppi_status=preparation['ppi_status'], best_r0=None,
        interpretation='Sensitivity experiment; no biological optimum selected. Compare identical LR panels, budgets, and spatial support.'))
    print(f'Complete: {out}/all_conditions.csv; no best r0 selected')


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fetch-lr'); f.add_argument('--output', default='data/cccvelo_lr')
    d = sub.add_parser('download-data'); d.add_argument('--output', default='data/mousebrain_bin60.h5ad')
    d.add_argument('--url', default=BRAIN_URL)
    i = sub.add_parser('inspect'); i.add_argument('--h5ad', required=True); i.add_argument('--output')
    a = sub.add_parser('prepare')
    a.add_argument('--h5ad', required=True); a.add_argument('--output', required=True)
    a.add_argument('--lr-dir', default='data/cccvelo_lr')
    a.add_argument('--matrix', required=True, help='X, raw, or layer:NAME; choose after inspect')
    a.add_argument('--expression-state', choices=['counts', 'processed'], required=True)
    a.add_argument('--spatial-key', default='spatial'); a.add_argument('--coordinate-unit', required=True)
    a.add_argument('--gene-key'); a.add_argument('--gene-map', help='One-to-one CSV: input_gene,lr_gene')
    a.add_argument('--symbol-policy', choices=['exact', 'author-uppercase'], default='exact')
    a.add_argument('--center', type=float, nargs=2); a.add_argument('--max-spots', type=int, default=128)
    a.add_argument('--pairs-per-type', type=int, default=4)
    a.add_argument('--min-fraction', type=float, default=0.05); a.add_argument('--min-mean', type=float, default=0.05)
    a.add_argument('--contact-radius-h', type=float, default=2.0)
    a.add_argument('--ppi'); a.add_argument('--score-provenance')
    r = sub.add_parser('run'); r.add_argument('--input', required=True); r.add_argument('--output', required=True)
    r.add_argument('--sensitivity', action='store_true', help='Nine eps/rho settings instead of one')
    return p


def main():
    args = parser().parse_args()
    if args.command == 'fetch-lr':
        fetch_lr(args.output)
    elif args.command == 'download-data':
        print(f'Downloading brain data to {args.output}...', flush=True)
        download(args.url, args.output, hdf5=True)
        write_json(str(args.output) + '.download.json', dict(url=args.url, sha256=digest(args.output),
            default_url_source=BRAIN_SOURCE,
            checksum_source='Locally calculated, not an author-published authenticity checksum'))
        print(f'Downloaded {args.output}; run inspect next')
    elif args.command == 'inspect':
        inspect_data(args.h5ad, args.output)
    elif args.command == 'prepare':
        prepare(args)
    else:
        run_stage(args)


if __name__ == '__main__':
    main()
