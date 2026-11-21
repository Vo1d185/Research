"""Experiment 2, expanded (reports/stage2_ppi_03/plan.md): data for the COMMOT-paper datasets. No OT, no S.

  fetch-zenodo        download archives of the COMMOT paper's Zenodo record 7272562 by name; verify the
                      Zenodo MD5 and record SHA256 (an existing file is verified, never overwritten)
  fetch-cellchat      CellChatDB v2 (jinworks/CellChat, the commit already pinned for mouse) for
                      --species human | mouse -> cellchat_inventory format (interaction/complex/gene_symbols.csv)
  fetch-db-human      scSeqComm human tables at the commit pinned in evaluate_targets.py
                      (TF_PPR_KEGG_human, TF_PPR_REACTOME_human, TF_TG_TRRUSTv2) -> evaluate_targets format
  prepare-a           A line (plan section 3): Stage 2 input + full-transcriptome gene export, same rules as
                      stage2.py prepare / export-genes. --dataset breast | seqfish
  prepare-b           B line (plan section 4): the paper protocol. LR list = the paper's own filtered list
                      (uns['commot-cellchat-info']['df_ligrec'] of the paper's output h5ad, COMMOT's CellChatDB v1,
                      secreted only); cutoff 1000 coordinate units; paper target-gene expression (HVG; raw
                      counts for the mouse brain, log1p otherwise) and the paper's received signal (step 0a).
                      --dataset post | breast | seqfish

Environment: commot-research (numpy, scipy, pandas, h5py, anndata, scanpy, rdata). All paths are arguments;
existing outputs are never overwritten.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import time
import zipfile

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.spatial import cKDTree

from research import check

# ------------------------------------------------------------------ pinned sources
ZENODO_URL = 'https://zenodo.org/api/records/7272562/files/{name}/content'
ZENODO_MD5 = {   # Zenodo API "checksum" of record 7272562 (version 0.0.2), read 2026-11-21
    'downstream_gene_correlation.zip': '2945fc722a492ac048d80ac032ae399a',
    'visium-mouse-brain.zip': '17b161cfa13f3e5af189eb32d05adaec',
    'visium-human-breast_cancer.zip': 'a48416fd41bd9a2bfeddbeedde452833',
    'seqfishplus-mouse-sscortex.zip': 'e610e9fdf954b3351058b9352944249a',
}
BREAST_SAMPLE = 'V1_Breast_Cancer_Block_A_Section_1'
BREAST_MEMBERS = {   # 10x files inside downstream_gene_correlation.zip (SHA256 computed locally 2026-11-21)
    'filtered_feature_bc_matrix.h5': (
        f'downstream_gene_correlation/data/{BREAST_SAMPLE}/filtered_feature_bc_matrix.h5',
        'b3cbb6ead60afad848ebcc07a34835613d9d2d79e416ad8300b210d9de584b65'),
    'spatial.tar.gz': (
        f'downstream_gene_correlation/data/{BREAST_SAMPLE}/{BREAST_SAMPLE}_spatial.tar.gz',
        '415530a9d31172f63fd358d5d1cefb0798698ef9c0e614f52825894cd45c61bf'),
}
SEQFISH_MEMBERS = {  # inside seqfishplus-mouse-sscortex.zip: the 523 cortex cells used by the paper
    'expression.npz': ('seqfishplus-mouse-sscortex/data/original_data/expression.npz',
                       '685bd690176f36baa80a1f978917719b16ffc4f246d28bbfc3384ed69f93ae37'),
    'genes.txt': ('seqfishplus-mouse-sscortex/data/original_data/genes.txt',
                  '732bf79e4c353a3054e9cf3159e3f5f79ca4f8b11c02e8ba293e30c1638bd4d5'),
    'positions.csv': ('seqfishplus-mouse-sscortex/data/original_data/positions.csv',
                      '9f538c7362a5decaeb7609e254dadc2938c97ae73cced384b839fcb7284d432f'),
}
SEQFISH_UM_PER_PX = 0.103      # Eng et al., Nature 568, 235 (2019): "X-Y coordinates are in pixels (103 nm per pixel)"
CONTACT_NN_MULTIPLE = 2.0      # seqFISH+ contact pairs: d <= 2 x median nearest-neighbour distance (plan 3.2)
PAPER_OUTPUT = {               # the paper's COMMOT output (dis_thr = 1000) per dataset
    'post': ('visium-mouse-brain.zip', 'visium-mouse-brain/data/adata_filteredcellchat_dis1000.h5ad'),
    'breast': ('visium-human-breast_cancer.zip', 'visium-human-breast_cancer/data/adata_filteredcellchat_dis1000.h5ad'),
    'seqfish': ('seqfishplus-mouse-sscortex.zip', 'seqfishplus-mouse-sscortex/data/adata_dis1000.h5ad'),
}
PAPER = {   # Zenodo notebooks 1-lr_signaling.ipynb / cell-level_downstream_gene_correlation.ipynb; paper text (medians)
    'post': dict(species='mouse', dis_thr=1000.0, min_cell_pct=0.05, hvg_min_disp=0.5, eval_scale='raw', input_scale='normalize_log1p',
                 n_pairs=250, n_eval=164, n_hvg=2793, mean=0.19681672570590536, median=0.180),
    'breast': dict(species='human', dis_thr=1000.0, min_cell_pct=0.05, hvg_min_disp=0.5, eval_scale='log1p', input_scale='normalize_log1p',
                   n_pairs=300, n_eval=192, n_hvg=3068, mean=0.23194032950633572, median=0.237),
    'seqfish': dict(species='mouse', dis_thr=1000.0, min_cell_pct=0.0, hvg_min_disp=0.1, eval_scale='log1p', input_scale='log1p',
                    n_pairs=598, n_eval=519, n_hvg=4797, mean=0.2312192487386623, median=0.230),
}
# input_scale: expression given to COMMOT in the notebook that wrote the dis1000 file. seqFISH+: cell 12 of
# seqfishplus-mouse-sscortex/1-lr_signaling.ipynb re-reads the counts and applies log1p only (no normalize_total).
PAPER_X_TOL = 1e-4             # prepare-b: our LR-gene expression vs X of the paper's output (float32)
HVG = dict(min_cells=3, min_mean=0.0125, max_mean=3.0)
CELLCHAT_SHA256 = {'mouse': 'dd59b19646e546b3ed41f394ad2cf64cf44f3d3fcfbb58ac176f0fe1322f9413',
                   'human': '582e99db2bb1b1e05c0364553893d804952bd007fef870720a9cc659a7bd3ad8'}
SCSEQCOMM_HUMAN = {   # scSeqComm 95c9bb2 (evaluate_targets.SCSEQCOMM_COMMIT); SHA256 computed locally 2026-11-21
    'TF_TG_TRRUSTv2.RData': '24b7667c6cce567aa6f3b19c0d84138f7a392a4784a409c72816d4fa43ea4445',
    'TF_PPR_KEGG_human.rda': 'ad0b72481ea50c3973a119e3d210342e9def7ece4159508dc0102666b5dc168a',
    'TF_PPR_REACTOME_human.rda': '1e935075050e2bc9479e79463c5e11aa8ba39baf38c57343183b8b5b37df4e60',
}
A_CONFIG = dict(eps=0.2, rho=1.0, r0_um=[50.0], diffusion_cutoff_um=363.0757, nitermax=50000, stopthr=1e-11)
B_CONFIG = dict(eps=0.1, rho=10.0, r0_um=[50.0], nitermax=50000, stopthr=1e-11)


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as s:
        for block in iter(lambda: s.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def md5(path):
    h = hashlib.md5()
    with Path(path).open('rb') as s:
        for block in iter(lambda: s.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False,
                                     default=lambda x: x.item() if hasattr(x, 'item') else str(x)))


def fresh(path):
    path = Path(path)
    check(not path.exists(), f'Output exists: {path}')
    path.mkdir(parents=True)
    return path


# ------------------------------------------------------------------ downloads
def fetch_zenodo(args):
    from stage1 import download
    folder = Path(args.output)
    folder.mkdir(parents=True, exist_ok=True)
    for name in args.name:
        check(name in ZENODO_MD5, f'Unknown archive {name}')
        target = folder/name
        if not target.exists():
            for attempt in range(3):
                try:
                    download(ZENODO_URL.format(name=name), target)
                    break
                except Exception as error:   # transient HTTP errors; download() leaves no partial file
                    if attempt == 2:
                        raise
                    print(f'retry {name}: {error}', flush=True)
                    time.sleep(10 * 2**attempt)
        check(md5(target) == ZENODO_MD5[name], f'{target}: MD5 differs from the Zenodo record')
        record = folder/(name + '.json')
        info = dict(name=name, md5=ZENODO_MD5[name], sha256=digest(target), size=target.stat().st_size,
                    url=ZENODO_URL.format(name=name), checked_utc=now())
        if not record.exists():
            write_json(record, info)
        print(json.dumps(info), flush=True)


def fetch_cellchat(args):
    """cellchat_inventory.fetch for either species (same commit, same conversion)."""
    import rdata
    from stage1 import download
    from cellchat_inventory import CELLCHAT_COMMIT
    out = Path(args.output)
    check(not (out/'db.json').exists(), f'{out} already has a converted database')
    url = f'https://raw.githubusercontent.com/jinworks/CellChat/{CELLCHAT_COMMIT}/data/CellChatDB.{args.species}.rda'
    raw = out/'raw'/f'CellChatDB.{args.species}.rda'
    download(url, raw, expected=CELLCHAT_SHA256[args.species])
    obj = next(iter(rdata.read_rda(raw).values()))
    interaction = obj['interaction'].copy()
    interaction.columns = [str(c) for c in interaction.columns]
    complex_ = obj['complex'].copy()
    complex_.index = complex_.index.astype(str)
    complex_.columns = [str(c) for c in complex_.columns]
    genes = obj['geneInfo'][['Symbol']].astype(str)
    interaction.to_csv(out/'interaction.csv', index=False)
    complex_.to_csv(out/'complex.csv', index_label='complex')
    genes.to_csv(out/'gene_symbols.csv', index=False)
    write_json(out/'db.json', dict(status='complete', source='CellChat v2 (jinworks/CellChat)', species=args.species,
        commit=CELLCHAT_COMMIT, url=url, raw_sha256=CELLCHAT_SHA256[args.species],
        sha256_source='calculated locally, not author-published', license='GPL-3 (downloaded, not redistributed)',
        retrieved_utc=now(), n_interactions=len(interaction),
        annotation_counts=interaction.annotation.value_counts().to_dict(),
        converted_sha256={f: digest(out/f) for f in ['interaction.csv', 'complex.csv', 'gene_symbols.csv']}))
    print(f'CellChatDB {args.species}: {len(interaction)} interactions -> {out}')


def load_cellchat(folder):
    """cellchat_inventory.load_cellchat for either species (the mouse file is the one already pinned there)."""
    folder = Path(folder)
    meta = json.loads((folder/'db.json').read_text())
    species = meta.get('species', 'mouse')
    check(meta.get('status') == 'complete' and meta['raw_sha256'] == CELLCHAT_SHA256[species], 'Run fetch-cellchat first')
    for name, sha in meta['converted_sha256'].items():
        check(digest(folder/name) == sha, f'{name} changed after conversion')
    it = pd.read_csv(folder/'interaction.csv', dtype=str)
    complexes = set(pd.read_csv(folder/'complex.csv', dtype=str).complex)
    symbols = set(pd.read_csv(folder/'gene_symbols.csv', dtype=str).Symbol)
    for side in ['ligand', 'receptor']:
        unknown = ~it[side].isin(complexes) & ~it[side].isin(symbols)
        check(not unknown.any(), f'{side} entities neither complex nor gene: {it.loc[unknown, side].head().tolist()}')
    check(not it.duplicated(['ligand', 'receptor']).any(), 'Duplicate CellChat ligand-receptor pairs')
    return it, dict(meta, species=species)


def cellchat_pairs(folder):
    """stage2.cellchat_pairs with the species-aware loader."""
    import stage2
    it, meta = load_cellchat(folder)
    cx = pd.read_csv(Path(folder)/'complex.csv', dtype=str)
    cols = [c for c in cx.columns if c.startswith('subunit')]
    complexes = {r.complex: [s for s in (getattr(r, c) for c in cols) if isinstance(s, str) and s]
                 for r in cx.itertuples(index=False)}
    pairs = it[it.annotation.isin(stage2.CLASSES)].copy()
    pairs['type'] = pairs.annotation.map(stage2.CLASSES)
    subunits = {e: complexes.get(e, [e]) for e in set(pairs.ligand) | set(pairs.receptor)}
    return pairs, subunits, meta


def fetch_db_human(args):
    """evaluate_targets.fetch_db for the human scSeqComm tables (KEGG + REACTOME PPR; TF_TG_TRRUSTv2, the
    library COMMOT used for human)."""
    import rdata
    from urllib.error import HTTPError
    from stage1 import download
    from evaluate_targets import SCSEQCOMM_COMMIT, SCSEQCOMM_RAW
    out = Path(args.output)
    check(not (out/'db.json').exists(), f'{out} already has a converted database')
    for name, sha in SCSEQCOMM_HUMAN.items():
        for attempt in range(4):
            try:
                download(SCSEQCOMM_RAW + name, out/'raw'/name, expected=sha)
                break
            except HTTPError as error:
                if error.code not in (429, 500, 502, 503, 504) or attempt == 3:
                    raise
                time.sleep(5 * 2**attempt)
    tf_tg = next(iter(rdata.read_rda(out/'raw'/'TF_TG_TRRUSTv2.RData').values()))
    edges = pd.DataFrame([(str(tf), str(t)) for tf, targets in tf_tg.items() for t in np.atleast_1d(targets)],
                         columns=['tf', 'target']).drop_duplicates()
    frames = []
    for name, source in [('TF_PPR_KEGG_human.rda', 'KEGG'), ('TF_PPR_REACTOME_human.rda', 'REACTOME')]:
        frame = next(iter(rdata.read_rda(out/'raw'/name).values()))
        frames.append(pd.DataFrame(dict(receptor=frame.receptor.astype(str), pathway=frame.pathway.astype(str),
                                        tf=frame.tf.astype(str), weight=frame.tf_PPR.astype(float), source=source)))
    receptor_tf = pd.concat(frames, ignore_index=True)
    edges.to_csv(out/'tf_tg.csv', index=False)
    receptor_tf.to_csv(out/'receptor_tf.csv', index=False)
    write_json(out/'db.json', dict(status='complete', source=f'scSeqComm {SCSEQCOMM_COMMIT}', species='human',
        license='GPL-3 (downloaded, not redistributed)', raw_sha256=SCSEQCOMM_HUMAN,
        weight='personalized PageRank of TF seeded at the receptor (scSeqComm TF_PPR)',
        tf_tg_edges=len(edges), tfs=int(edges.tf.nunique()), receptor_tf_rows=len(receptor_tf),
        receptors=int(receptor_tf.receptor.nunique()),
        converted_sha256={f: digest(out/f) for f in ['tf_tg.csv', 'receptor_tf.csv']}))
    print(f'Converted scSeqComm human tables -> {out}')


# ------------------------------------------------------------------ raw data readers
def zip_member(zip_path, member, sha):
    with zipfile.ZipFile(zip_path) as z:
        data = z.read(member)
    check(hashlib.sha256(data).hexdigest() == sha, f'{zip_path}:{member} checksum mismatch')
    return data


def register_breast():
    import stage2
    files = {k: sha for k, (_, sha) in BREAST_MEMBERS.items()}
    stage2.SAMPLES.setdefault('breast', (BREAST_SAMPLE, files))
    check(stage2.SAMPLES['breast'] == (BREAST_SAMPLE, files), 'stage2.SAMPLES["breast"] differs')


def breast_folder(zenodo_dir, folder):
    """Write the two 10x files of the breast-cancer slide (from the paper's Zenodo copy) under 10x names."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    for suffix, (member, sha) in BREAST_MEMBERS.items():
        target = folder/f'{BREAST_SAMPLE}_{suffix}'
        if target.exists():
            check(digest(target) == sha, f'{target} exists with other content')
        else:
            target.write_bytes(zip_member(Path(zenodo_dir)/'downstream_gene_correlation.zip', member, sha))
    return folder


def read_dataset(dataset, zenodo_dir, data_dir):
    """Counts (cells/spots x genes, CSR float64), unique var names, ids, micrometre positions, raw coordinates,
    micrometres per coordinate unit, nearest-neighbour spacing (um), array indices (Visium) and provenance."""
    import stage2
    if dataset in ('post', 'breast'):
        key = 'posterior' if dataset == 'post' else 'breast'
        if dataset == 'breast':
            register_breast()
            breast_folder(zenodo_dir, data_dir)
        v = stage2.read_visium(data_dir, key)
        xy, um_per_px, pitch = stage2.micrometres(v['pos'], v['scalefactors'])
        raw_xy = v['pos'][['pxl_col', 'pxl_row']].to_numpy(float)
        return dict(X=v['X'], var=v['var'], ids=pd.Index(v['barcodes'].astype(str)), xy=xy, raw_xy=raw_xy,
                    um_per_unit=um_per_px, spacing_um=pitch, array=v['pos'][['array_row', 'array_col']].to_numpy(),
                    sample=v['sample'], unit='full-resolution pixel', n_duplicated=v['n_duplicated'],
                    source=dict(folder=str(Path(data_dir).resolve()),
                                files={s: sha for s, sha in stage2.SAMPLES[key][1].items()}))
    check(dataset == 'seqfish', f'Unknown dataset {dataset}')
    from anndata.utils import make_index_unique
    z = Path(zenodo_dir)/'seqfishplus-mouse-sscortex.zip'
    blobs = {k: zip_member(z, m, sha) for k, (m, sha) in SEQFISH_MEMBERS.items()}
    X = sparse.load_npz(io.BytesIO(blobs['expression.npz'])).tocsr().astype(np.float64)
    genes = np.loadtxt(io.BytesIO(blobs['genes.txt']), dtype=str)
    pts = np.loadtxt(io.BytesIO(blobs['positions.csv']))
    check(X.shape == (len(pts), len(genes)) == (523, 10000), f'Unexpected seqFISH+ shapes {X.shape}')
    check(np.allclose(X.data, np.round(X.data)) and X.data.min() >= 0, 'seqFISH+ expression must be counts')
    var = pd.Index(make_index_unique(pd.Index(genes)))
    xy = pts * SEQFISH_UM_PER_PX
    nn = cKDTree(xy).query(xy, k=2)[0][:, 1]
    check((nn > 0).all(), 'Coincident cells')
    return dict(X=X, var=var, ids=pd.Index([f'cell_{i:03d}' for i in range(len(pts))]), xy=xy, raw_xy=pts,
                um_per_unit=SEQFISH_UM_PER_PX, spacing_um=float(np.median(nn)), array=None,
                sample='seqFISH+ mouse SS cortex (Eng et al. 2019; COMMOT Zenodo original_data)', unit='pixel',
                n_duplicated=int(pd.Index(genes).duplicated().sum()),
                source=dict(archive=str(z.resolve()), archive_md5=ZENODO_MD5[z.name],
                            members={k: sha for k, (_, sha) in SEQFISH_MEMBERS.items()}))


def distance_pairs(xy, cutoff):
    """Ordered pairs (i, j), i != j, with d(i, j) <= cutoff."""
    pairs = cKDTree(xy).query_pairs(cutoff, output_type='ndarray')
    s = np.concatenate([pairs[:, 0], pairs[:, 1]]).astype(np.int64)
    t = np.concatenate([pairs[:, 1], pairs[:, 0]]).astype(np.int64)
    keep = np.linalg.norm(xy[s] - xy[t], axis=1) <= cutoff
    return s[keep], t[keep]


def entity_expression(Y, var, entities, subunits):
    import stage2
    return stage2.entity_expression(Y, var, entities, subunits)


# ------------------------------------------------------------------ gene exports
def gene_export(out, X_counts, var, ids, inp, hvg_min_disp=0.5):
    """stage2.export_genes for a matrix in memory: normalize_total + log1p of every gene, entities re-derived and
    compared with expression.csv, COMMOT HVG flag; files gene_expression.npz, genes.csv, gene_export.json."""
    import stage2
    out = fresh(out)
    Y, target = stage2.normalize_log1p(X_counts)
    X = Y.toarray()
    expression = pd.read_csv(Path(inp)/'expression.csv', index_col=0)
    check(list(expression.index.astype(str)) == list(ids), 'Gene export ids differ from expression.csv')
    ent = pd.read_csv(Path(inp)/'entities.csv', dtype=str)
    pos = {g: i for i, g in enumerate(var)}
    for r in ent.itertuples(index=False):
        derived = X[:, [pos[g] for g in r.subunits.split(';')]].min(axis=1)
        check(np.allclose(derived, expression[r.entity].to_numpy(float), rtol=1e-9, atol=1e-12),
              f'{r.entity} differs from expression.csv')
    hvg = hvg_flags(X_counts, var, hvg_min_disp)
    np.savez_compressed(out/'gene_expression.npz', X=X.astype(np.float32), spots=np.asarray(ids, str),
                        genes=np.asarray(var, dtype=str))
    pd.DataFrame(dict(data_symbol=var, lr_symbol=var, detected_fraction=(X > 0).mean(axis=0),
                      commot_hvg=hvg.to_numpy())).to_csv(out/'genes.csv', index=False)
    write_json(out/'gene_export.json', dict(status='complete', created_utc=now(), input=str(Path(inp).resolve()),
        normalization=f'normalize_total (target = median total count = {target:g}) then log1p',
        commot_hvg=dict(HVG, min_disp=hvg_min_disp, n=int(hvg.sum()), flavor='seurat (scanpy default)'),
        gene_expression_sha256=digest(out/'gene_expression.npz'), spots=len(ids), genes=len(var),
        note='lr_symbol = data symbol (CellChat symbols of the same species)'))
    return out


def hvg_flags(X_counts, var, min_disp):
    """COMMOT's downstream-gene notebooks: filter_genes(min_cells=3) on counts, log1p, highly_variable_genes
    (min_mean=0.0125, max_mean=3, min_disp; seurat flavour). Genes removed by the filter are not HVG."""
    import anndata
    import scanpy as sc
    adata = anndata.AnnData(X=sparse.csr_matrix(X_counts, dtype=np.float32), var=pd.DataFrame(index=pd.Index(var)))
    sc.pp.filter_genes(adata, min_cells=HVG['min_cells'])
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, min_mean=HVG['min_mean'], max_mean=HVG['max_mean'], min_disp=min_disp)
    return pd.Series(adata.var.highly_variable.reindex(pd.Index(var)).fillna(False).astype(bool).to_numpy(),
                     index=pd.Index(var))


# ------------------------------------------------------------------ A line
def prepare_a(args):
    import stage2
    out = Path(args.output)
    check(not out.exists() and not Path(args.genes_output).exists(), 'Output exists')
    d = read_dataset(args.dataset, args.zenodo_dir, args.data)
    Y, target = stage2.normalize_log1p(d['X'])
    frac = np.asarray((d['X'] > 0).sum(axis=0)).ravel() / d['X'].shape[0]
    pairs, subunits, cmeta = cellchat_pairs(args.cellchat)
    species = 'human' if args.dataset == 'breast' else 'mouse'
    check(cmeta['species'] == species, f'{args.dataset} needs the {species} CellChatDB, got {cmeta["species"]}')
    audit = stage2.commot_filter(pairs, subunits, d['var'], frac, args.min_cell_pct)
    selected = audit[audit.reason == 'selected'].copy()
    check(len(selected) and set(selected.type) == {'diffusion', 'contact'}, 'Both LR types must remain')
    selected = selected.sort_values(['type', 'ligand', 'receptor']).reset_index(drop=True)
    entities = sorted(set(selected.ligand) | set(selected.receptor))
    expr = entity_expression(Y, d['var'], entities, subunits)
    expr.index = d['ids']
    spatial = pd.DataFrame(dict(x=d['xy'][:, 0], y=d['xy'][:, 1]), index=d['ids'])
    if d['array'] is not None:
        spatial['array_row'], spatial['array_col'] = d['array'][:, 0], d['array'][:, 1]
        contact = dict(contact_rule='hex', contact='self + hexagonal neighbours (array_row/array_col)')
    else:
        cut = CONTACT_NN_MULTIPLE * d['spacing_um']
        contact = dict(contact_rule='distance', contact_cutoff_um=cut,
                       contact=f'self + cells within {CONTACT_NN_MULTIPLE:g} x median nearest-neighbour distance '
                               f'({d["spacing_um"]:.3f} um) = {cut:.3f} um')
    out.mkdir(parents=True)
    expr.to_csv(out/'expression.csv', index_label='spot')
    spatial.to_csv(out/'spatial.csv', index_label='spot')
    selected[['ligand', 'receptor', 'type', 'annotation', 'pathway_name', 'interaction_name']].to_csv(out/'lr.csv', index=False)
    pd.DataFrame(dict(entity=entities, subunits=[';'.join(subunits[e]) for e in entities])).to_csv(out/'entities.csv', index=False)
    audit.drop(columns=['ligand_complex', 'receptor_complex'], errors='ignore')[
        ['ligand', 'receptor', 'type', 'pathway_name', 'reason', 'min_subunit_fraction']].to_csv(out/'filter_audit.csv', index=False)
    cfg = dict(A_CONFIG, min_cell_pct=args.min_cell_pct, coordinate_unit='micrometer',
               expression_level='gene; complexes = subunit minimum', expression_state='processed', allow_self=True,
               dataset=args.dataset, species=species, line='A', plan='reports/stage2_ppi_03/plan.md', **contact)
    write_json(out/'config.json', cfg)
    n_contact = None
    if contact['contact_rule'] == 'distance':
        s, _ = distance_pairs(d['xy'], contact['contact_cutoff_um'])
        n_contact = np.bincount(s, minlength=len(d['ids']))
    files = ['expression.csv', 'spatial.csv', 'lr.csv', 'entities.csv', 'config.json']
    write_json(out/'preparation.json', dict(status='complete', created_utc=now(), dataset=args.dataset, line='A',
        sample=d['sample'], species=species, source=d['source'], cellchat=dict(species=cmeta['species'],
        commit=cmeta['commit'], raw_sha256=cmeta['raw_sha256']),
        normalization=f'normalize_total (target = median total count = {target:g}) then log1p',
        filter=f'COMMOT filter_lr_database rule: each subunit detected in >= {args.min_cell_pct:g} of all {d["X"].shape[0]} spots/cells',
        coordinate_unit=d['unit'], um_per_unit=d['um_per_unit'], neighbour_spacing_um=d['spacing_um'],
        contact_neighbours_per_cell=None if n_contact is None else dict(
            median=float(np.median(n_contact)), min=int(n_contact.min()), max=int(n_contact.max()),
            zero=int((n_contact == 0).sum())),
        spots=len(d['ids']), genes=len(d['var']), duplicated_features=d['n_duplicated'],
        pairs_selected={k: int(v) for k, v in selected.type.value_counts().items()},
        ligands=int(selected.ligand.nunique()), receptors=int(selected.receptor.nunique()),
        filter_counts={f'{t}/{r}': int(c) for (t, r), c in audit.groupby(['type', 'reason']).size().items()},
        prepared_sha256={f: digest(out/f) for f in files}))
    gene_export(args.genes_output, d['X'], d['var'], d['ids'], out)
    print(f'{args.dataset} A: {len(d["ids"])} spots, {len(selected)} pairs {dict(selected.type.value_counts())} -> {out}; '
          f'genes -> {args.genes_output}', flush=True)


# ------------------------------------------------------------------ B line (paper protocol)
def paper_output(dataset, zenodo_dir, work):
    """Path of the paper's output h5ad, extracted once from its Zenodo archive into --work."""
    name, member = PAPER_OUTPUT[dataset]
    target = Path(work)/Path(member).name.replace('.h5ad', f'_{dataset}.h5ad')
    if not target.exists():
        Path(work).mkdir(parents=True, exist_ok=True)
        archive = Path(zenodo_dir)/name
        check(archive.is_file() and md5(archive) == ZENODO_MD5[name], f'{archive} missing or not the Zenodo file')
        part = target.with_name(target.name + '.part')
        with zipfile.ZipFile(archive) as z, z.open(member) as src, part.open('wb') as dst:
            for block in iter(lambda: src.read(1 << 22), b''):
                dst.write(block)
        part.rename(target)
    return target


def prepare_b(args):
    import anndata
    import stage2
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    p = PAPER[args.dataset]
    d = read_dataset(args.dataset, args.zenodo_dir, args.data)
    h5 = paper_output(args.dataset, args.zenodo_dir, args.work)
    a = anndata.read_h5ad(h5)
    info = a.uns['commot-cellchat-info']
    df = info['df_ligrec']
    df = df if isinstance(df, pd.DataFrame) else pd.DataFrame(df)
    df.columns = [str(c) for c in df.columns]
    check(list(df.columns[:2]) == ['ligand', 'receptor'], f'Unexpected df_ligrec columns {list(df.columns)}')
    check(float(np.asarray(info['distance_threshold']).ravel()[0]) == p['dis_thr'], 'Paper distance threshold differs')
    check(len(df) == p['n_pairs'] and not df.duplicated(['ligand', 'receptor']).any(),
          f'{args.dataset}: paper LR list has {len(df)} pairs, expected {p["n_pairs"]}')
    # spot order: the paper's obs vs ours
    obs = pd.Index(a.obs_names.astype(str))
    if args.dataset == 'seqfish':
        check(list(obs) == [str(i) for i in range(len(d['ids']))], 'Paper seqFISH+ obs are not 0..n-1 in order')
        order = np.arange(len(obs))
    else:
        order = obs.get_indexer(d['ids'])
        check((order >= 0).all() and len(obs) == len(d['ids']), 'Paper barcodes differ from the 10x matrix')
    paper_xy = np.asarray(a.obsm['spatial'], float)[order]
    check(np.allclose(paper_xy, d['raw_xy']) or np.allclose(np.abs(paper_xy), np.abs(d['raw_xy'])),
          'Paper coordinates differ from the raw coordinates')
    if p['input_scale'] == 'log1p':
        Y, target = d['X'].copy().tocsr(), None
        Y.data = np.log1p(Y.data)
    else:
        Y, target = stage2.normalize_log1p(d['X'])
    subunits = {e: e.split('_') for e in set(df.ligand) | set(df.receptor)}
    missing = sorted({g for s in subunits.values() for g in s if g not in set(d['var'])})
    check(not missing, f'Paper LR subunits not measured: {missing[:5]}')
    lr_genes = sorted({g for s in subunits.values() for g in s})
    mine = Y[:, pd.Index(d['var']).get_indexer(lr_genes)].toarray()
    pidx = pd.Index(a.var_names.astype(str)).get_indexer(lr_genes)
    check((pidx >= 0).all(), 'LR genes missing from the paper output var')
    theirs = a.X[order][:, pidx]
    theirs = theirs.toarray() if sparse.issparse(theirs) else np.asarray(theirs)
    x_diff = float(np.abs(mine - theirs).max())
    check(x_diff <= PAPER_X_TOL, f'Our expression differs from the paper input X (max |diff| {x_diff:.3g})')
    entities = sorted(subunits)
    expr = entity_expression(Y, d['var'], entities, subunits)
    expr.index = d['ids']
    # COMMOT's own entity expression is X with heteromeric 'min'; compare with the paper's input X where possible
    lr = pd.DataFrame(dict(ligand=df.ligand.astype(str), receptor=df.receptor.astype(str), type='diffusion',
                           annotation='Secreted Signaling',
                           pathway_name=df.iloc[:, 2].astype(str) if df.shape[1] > 2 else '',
                           interaction_name=df.ligand.astype(str) + '_' + df.receptor.astype(str)))
    spatial = pd.DataFrame(dict(x=d['xy'][:, 0], y=d['xy'][:, 1]), index=d['ids'])
    if d['array'] is not None:
        spatial['array_row'], spatial['array_col'] = d['array'][:, 0], d['array'][:, 1]
    cutoff_um = p['dis_thr'] * d['um_per_unit']
    out.mkdir(parents=True)
    expr.to_csv(out/'expression.csv', index_label='spot')
    spatial.to_csv(out/'spatial.csv', index_label='spot')
    lr.to_csv(out/'lr.csv', index=False)
    pd.DataFrame(dict(entity=entities, subunits=[';'.join(subunits[e]) for e in entities])).to_csv(out/'entities.csv', index=False)
    cfg = dict(B_CONFIG, diffusion_cutoff_um=cutoff_um, paper_dis_thr=p['dis_thr'], coordinate_unit='micrometer',
               expression_level='gene; complexes = subunit minimum (COMMOT heteromeric_rule min)',
               expression_state='processed', allow_self=True, contact_rule='none',
               contact='no contact pairs (COMMOT secreted pairs only)', dataset=args.dataset, species=p['species'],
               line='B', plan='reports/stage2_ppi_03/plan.md')
    write_json(out/'config.json', cfg)
    # paper target-gene expression (cell-level_downstream_gene_correlation.ipynb) and the paper's received signal
    hvg = hvg_flags(d['X'], d['var'], p['hvg_min_disp'])
    cols = np.flatnonzero(hvg.to_numpy())
    E = d['X'][:, cols].toarray()
    if p['eval_scale'] == 'log1p':
        E = np.log1p(E)
    rec = a.obsm['commot-cellchat-sum-receiver']
    rec = rec if isinstance(rec, pd.DataFrame) else pd.DataFrame(rec)
    names = [f'r-{l}-{r}' for l, r in zip(lr.ligand, lr.receptor)]
    check(set(names) <= set(rec.columns), 'Paper received signal lacks some LR pairs')
    signal = rec[names].to_numpy(float)[order].T
    np.savez_compressed(out/'paper_eval.npz', E=E.astype(np.float32), genes=np.asarray(d['var'][cols], str),
                        spots=np.asarray(d['ids'], str), paper_signal=signal, ligand=lr.ligand.to_numpy(str),
                        receptor=lr.receptor.to_numpy(str))
    files = ['expression.csv', 'spatial.csv', 'lr.csv', 'entities.csv', 'config.json']
    write_json(out/'preparation.json', dict(status='complete', created_utc=now(), dataset=args.dataset, line='B',
        sample=d['sample'], species=p['species'], source=d['source'], paper_output=dict(file=str(h5.resolve()),
        archive=PAPER_OUTPUT[args.dataset][0], member=PAPER_OUTPUT[args.dataset][1], sha256=digest(h5)),
        paper=p, normalization=(f'normalize_total (target = median total count = {target:g}) then log1p' if target
                                else 'log1p of counts (no normalize_total), as the paper notebook'),
        paper_input_x_max_abs_diff=x_diff, paper_input_x_genes=len(lr_genes),
        coordinate_unit=d['unit'], um_per_unit=d['um_per_unit'], neighbour_spacing_um=d['spacing_um'],
        diffusion_cutoff_um=cutoff_um, spots=len(d['ids']), genes=len(d['var']),
        pairs=len(lr), ligands=int(lr.ligand.nunique()), receptors=int(lr.receptor.nunique()),
        paper_eval=dict(scale=p['eval_scale'], hvg_min_disp=p['hvg_min_disp'], n_hvg=int(len(cols)),
                        n_hvg_paper=p['n_hvg'], n_hvg_matches=bool(len(cols) == p['n_hvg']),
                        sha256=digest(out/'paper_eval.npz')),
        prepared_sha256={f: digest(out/f) for f in files}))
    print(f'{args.dataset} B: {len(d["ids"])} spots, {len(lr)} paper pairs, cutoff {cutoff_um:.4f} um, '
          f'{len(cols)} HVG (paper {p["n_hvg"]}) -> {out}', flush=True)


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    z = sub.add_parser('fetch-zenodo'); z.add_argument('--name', nargs='+', required=True); z.add_argument('--output', required=True)
    c = sub.add_parser('fetch-cellchat'); c.add_argument('--species', choices=['human', 'mouse'], required=True)
    c.add_argument('--output', required=True)
    h = sub.add_parser('fetch-db-human'); h.add_argument('--output', required=True)
    a = sub.add_parser('prepare-a'); a.add_argument('--dataset', choices=['breast', 'seqfish'], required=True)
    a.add_argument('--zenodo-dir', required=True); a.add_argument('--data', help='folder for the breast 10x files')
    a.add_argument('--cellchat', required=True); a.add_argument('--min-cell-pct', type=float, default=0.05)
    a.add_argument('--output', required=True); a.add_argument('--genes-output', required=True)
    b = sub.add_parser('prepare-b'); b.add_argument('--dataset', choices=['post', 'breast', 'seqfish'], required=True)
    b.add_argument('--zenodo-dir', required=True)
    b.add_argument('--data', help='10x folder: posterior = stage2 fetch folder; breast = the prepare-a --data folder')
    b.add_argument('--work', required=True, help='where the paper output h5ad is extracted (git-ignored)')
    b.add_argument('--output', required=True)
    return p


if __name__ == '__main__':
    args = parser().parse_args()
    {'fetch-zenodo': fetch_zenodo, 'fetch-cellchat': fetch_cellchat, 'fetch-db-human': fetch_db_human,
     'prepare-a': prepare_a, 'prepare-b': prepare_b}[args.command](args)
