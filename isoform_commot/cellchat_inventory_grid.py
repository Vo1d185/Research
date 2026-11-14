"""Extended CellChatDB inventory: dataset x threshold x target library x complex rule.

Counts only — no panel selection, no OT.

  fetch-cccvelo-db  CCCvelo Databases/Database.rar at commit fac213c (SHA256 pinned),
                    extracted with bsdtar; RecTFDB -> receptor_tf.csv, TFTGDB -> tf_tg.csv
  fetch-visium      10x Visium V1_Mouse_Brain_Sagittal_Posterior (Space Ranger 1.1.0), the
                    sample COMMOT's notebooks load via scanpy.datasets.visium_sge
  report            the grid below, for Secreted Signaling and Cell-Cell Contact

Steps (as in cellchat_inventory.py):
  0  all CellChat pairs of the class
  a  single_gene_only: ligand and receptor are single genes; complex rules: all pairs
  b  ligand and receptor detected in >= threshold of spots
  c1 receptor has a receptor->TF link in the target library
  c  top-10 TFs (max weight over the receptor's rows) -> union of their targets, minus
     every gene (incl. subunits) of the pairs passing (b) in either class, kept if
     detected in >= threshold of spots; >= 5 targets
     (the threshold applies to targets too, as CCCvelo's pct.ct filters L, R and TG genes)

Complex rules:
  single_gene_only         inventory 01 rule (complexes excluded at a)
  complex_any_subunit_tf   complex expression = min over subunits (COMMOT heteromeric_rule='min');
                           c1 if ANY subunit has a TF link; TFs ranked over all subunits
  complex_all_subunits_tf  same, but c1 needs EVERY subunit to have a TF link
Detection of a complex (main): fraction of spots where the subunit minimum is > 0, i.e. all
subunits detected in the same spot. Comparison: 'per_subunit' = every subunit separately
>= threshold, which is what COMMOT's filter_lr_database actually does (heteromeric, 'min').

Target libraries: scSeqComm mouse (KEGG+Reactome PPR; mouse symbols) and CCCvelo
RecTFDB/TFTGDB (upper-case symbols; matched to mouse symbols by upper-casing, which is
symbol matching, not ortholog mapping).
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import shutil
import subprocess

import numpy as np
import pandas as pd
from scipy import sparse

from cellchat_inventory import sharing
from research import check

CCC_COMMIT = 'fac213c600e7929b5d14e7d95dbdaabd722f9bc9'
DB_URL = f'https://raw.githubusercontent.com/SunXQlab/CCCvelo/{CCC_COMMIT}/Databases/Database.rar'
DB_SHA256 = '0afe64a0f0c7c5db14fba2a8c50e031d534233556a03254c46a2c1e9a3c016df'
VISIUM_SAMPLE = 'V1_Mouse_Brain_Sagittal_Posterior'
VISIUM_PREFIX = f'https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/{VISIUM_SAMPLE}'
VISIUM_FILES = {
    f'{VISIUM_SAMPLE}_filtered_feature_bc_matrix.h5': 'df962209143860488bd679ff7a1b3bb14985208d0929fb814934ba7571c488b6',
    f'{VISIUM_SAMPLE}_spatial.tar.gz': 'b878d7eaca487ffc17224a870f86599abf9d45a99033c1ef65661ff119cef0d5',
}
THRESHOLDS = (0.01, 0.05)
TOP_TF, MIN_TARGETS = 10, 5
COMPLEX_RULES = ('single_gene_only', 'complex_any_subunit_tf', 'complex_all_subunits_tf')
CLASSES = ['Secreted Signaling', 'Cell-Cell Contact']


# ------------------------------------------------------------------ fetch
def fetch_cccvelo_db(output, lr_dir):
    from stage1 import digest, download, load_lr, write_json
    out = Path(output)
    check(not (out/'db.json').exists(), f'{out} already has a converted database')
    rar = out/'raw'/'Database.rar'
    download(DB_URL, rar, expected=DB_SHA256)
    tool = shutil.which('bsdtar')
    check(tool is not None, 'Need bsdtar (libarchive >= 3.4 reads RAR5) to extract Database.rar')
    listing = subprocess.run([tool, '-tf', str(rar)], capture_output=True, text=True, check=True).stdout.split('\n')
    listing = [n for n in listing if n]
    check(all(not n.startswith('/') and '..' not in Path(n).parts for n in listing), 'Unsafe paths in archive')
    ext = out/'extracted'
    check(not ext.exists(), f'{ext} exists; remove it to re-extract')
    ext.mkdir(parents=True)
    subprocess.run([tool, '-xf', str(rar), '-C', str(ext)], check=True)
    contents = {str(p.relative_to(ext)): dict(bytes=p.stat().st_size, sha256=digest(p))
                for p in sorted(ext.rglob('*')) if p.is_file()}
    base = ext/'Database'
    rec = pd.read_csv(base/'RecTFDB.csv', dtype={'source': str, 'target': str})
    tftg = pd.read_csv(base/'TFTGDB.csv', dtype={'source': str, 'target': str})
    lig = pd.read_csv(base/'LigRecDB.csv', dtype=str, encoding='utf-8-sig')
    check({'source', 'target', 'score'} <= set(rec) and {'source', 'target', 'score'} <= set(tftg),
          'Unexpected RecTFDB/TFTGDB columns')
    check(not rec.duplicated(['source', 'target']).any() and not tftg.duplicated(['source', 'target']).any(),
          'Duplicate rows in RecTFDB/TFTGDB')
    rec.rename(columns={'source': 'receptor', 'target': 'tf', 'score': 'weight'})[
        ['receptor', 'tf', 'weight']].to_csv(out/'receptor_tf.csv', index=False)
    tftg.rename(columns={'source': 'tf', 'target': 'target', 'score': 'n_sources'})[
        ['tf', 'target', 'n_sources']].to_csv(out/'tf_tg.csv', index=False)
    pinned = load_lr(lr_dir)
    same_lr = set(zip(lig.source, lig.target)) == set(zip(pinned.ligand, pinned.receptor))
    write_json(out/'db.json', dict(status='complete', source=f'CCCvelo {CCC_COMMIT} Databases/Database.rar',
        url=DB_URL, raw_sha256=DB_SHA256, sha256_source='calculated locally at first download',
        extracted_with=subprocess.run([tool, '--version'], capture_output=True, text=True).stdout.strip(),
        contents=contents, receptor_tf_rows=len(rec), receptors=int(rec.source.nunique()),
        tfs=int(rec.target.nunique()), tf_tg_edges=len(tftg), tf_tg_tfs=int(tftg.source.nunique()),
        tf_tg_targets=int(tftg.target.nunique()),
        ligrecdb_equals_pinned_diff_plus_cont_tables=bool(same_lr),
        weight='RecTFDB score (receptor->TF); TFTGDB score kept as n_sources, not used',
        symbols='upper-case; matched to mouse symbols by upper-casing (not ortholog mapping)',
        retrieved_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'),
        converted_sha256={f: digest(out/f) for f in ['receptor_tf.csv', 'tf_tg.csv']}))
    print(f'CCCvelo database: {len(rec)} receptor-TF rows, {len(tftg)} TF-TG edges -> {out}')


def fetch_visium(output):
    from stage1 import digest, download, write_json
    out = Path(output)
    check(not (out/'source.json').exists(), f'{out} already fetched')
    for name, sha in VISIUM_FILES.items():
        download(f'{VISIUM_PREFIX}/{name}', out/name, expected=sha)
    with (out/f'{VISIUM_SAMPLE}_filtered_feature_bc_matrix.h5').open('rb') as stream:
        check(stream.read(8) == b'\x89HDF\r\n\x1a\n', 'Matrix file is not HDF5')
    write_json(out/'source.json', dict(status='complete', sample=VISIUM_SAMPLE, url_prefix=VISIUM_PREFIX,
        files={n: dict(sha256=s, bytes=(out/n).stat().st_size) for n, s in VISIUM_FILES.items()},
        sha256_source='calculated locally at first download, not published by 10x',
        why='COMMOT docs/notebooks (Basic_usage, visium-mouse_brain) at d117445 load this sample with '
            'scanpy.datasets.visium_sge, which uses Space Ranger 1.1.0 for V1_ samples',
        retrieved_utc=datetime.now(timezone.utc).isoformat(timespec='seconds')))
    print(f'Visium {VISIUM_SAMPLE} -> {out}')


# ------------------------------------------------------------------ data
def load_stereo(h5ad, matrix):
    import anndata as ad
    data = ad.read_h5ad(h5ad)
    X = data.layers[matrix[6:]] if matrix.startswith('layer:') else (data.raw.X if matrix == 'raw' else data.X)
    genes = pd.Index(data.var_names.astype(str))
    check(genes.is_unique, 'Duplicate gene symbols in h5ad')
    return dict(X=sparse.csc_matrix(X), genes=genes, n_spots=data.n_obs, duplicated_symbols=[])


def load_visium(h5):
    """10x filtered matrix -> spots x unique symbols; duplicated symbols are summed."""
    import h5py
    with h5py.File(h5, 'r') as f:
        g = f['matrix']
        M = sparse.csc_matrix((g['data'][:], g['indices'][:], g['indptr'][:]), shape=tuple(g['shape'][:]))
        names = np.array([x.decode() for x in g['features']['name'][:]])
        types = {x.decode() for x in g['features']['feature_type'][:]}
    check(types == {'Gene Expression'}, f'Unexpected feature types {types}')
    uniq, inverse = np.unique(names, return_inverse=True)
    A = sparse.csr_matrix((np.ones(len(names)), (np.arange(len(names)), inverse)), shape=(len(names), len(uniq)))
    X = (M.T.tocsr() @ A).tocsc()
    duplicated = sorted(set(names[pd.Series(names).duplicated().to_numpy()]))
    return dict(X=X, genes=pd.Index(uniq), n_spots=X.shape[0], duplicated_symbols=duplicated)


def gene_fractions(X, n_spots):
    return np.asarray((X > 0).sum(axis=0)).ravel() / n_spots


def entity_detection(subunits, X, genes, frac):
    """entity -> co-detection (all subunits > 0 in the same spot) and per-subunit minimum fraction."""
    pos = {g: i for i, g in enumerate(genes)}
    rows = {}
    for entity, subs in subunits.items():
        missing = [s for s in subs if s not in pos]
        if missing:
            rows[entity] = dict(codetect=np.nan, per_subunit=np.nan, missing=';'.join(missing))
            continue
        idx = [pos[s] for s in subs]
        if len(idx) == 1:
            value = float(frac[idx[0]])
            rows[entity] = dict(codetect=value, per_subunit=value, missing='')
            continue
        hits = np.asarray((X[:, idx] > 0).sum(axis=1)).ravel()
        rows[entity] = dict(codetect=float((hits == len(idx)).mean()),
                            per_subunit=float(min(frac[i] for i in idx)), missing='')
    return pd.DataFrame.from_dict(rows, orient='index')


# ------------------------------------------------------------------ target libraries
class Library:
    """receptor->TF (weighted) and TF->target tables in their own symbol space."""

    def __init__(self, name, receptor_tf, tf_tg, upper):
        self.name, self.upper = name, upper
        key = (lambda s: s.upper()) if upper else (lambda s: s)
        self.key = key
        rt = receptor_tf[['receptor', 'tf', 'weight']].copy()
        rt['receptor'], rt['tf'] = rt.receptor.map(key), rt.tf.map(key)
        self.by_receptor = {r: g[['tf', 'weight']] for r, g in rt.groupby('receptor')}
        tt = tf_tg[['tf', 'target']].assign(tf=tf_tg.tf.map(key), target=tf_tg.target.map(key))
        self.targets_of = tt.groupby('tf').target.apply(set).to_dict()
        self._cache = {}

    def receptor(self, subs):
        """(subunits with a link, top-TF list) for a receptor entity given as subunit symbols."""
        cache_key = tuple(subs)
        if cache_key not in self._cache:
            keys = [self.key(s) for s in subs]
            linked = [s for s, k in zip(subs, keys) if k in self.by_receptor]
            frames = [self.by_receptor[k] for k in keys if k in self.by_receptor]
            tfs = []
            if frames:
                ranked = pd.concat(frames).groupby('tf').weight.max().reset_index().sort_values(
                    ['weight', 'tf'], ascending=[False, True])
                tfs = ranked.tf.head(TOP_TF).tolist()
            self._cache[cache_key] = (linked, tfs)
        return self._cache[cache_key]

    def target_symbols(self, tfs, genes_by_key):
        keys = set().union(*[self.targets_of.get(t, set()) for t in tfs]) if tfs else set()
        return {genes_by_key[k] for k in keys if k in genes_by_key}


def genes_by_key(genes, upper):
    if not upper:
        return {g: g for g in genes}
    table = pd.Series(list(genes), index=[g.upper() for g in genes])
    return table[~table.index.duplicated(keep=False)].to_dict()  # drop case-collisions


# ------------------------------------------------------------------ grid
def subunit_map(it, complex_table):
    sub_cols = [c for c in complex_table.columns if c.startswith('subunit')]
    cx = {row.complex: [s for s in (getattr(row, c) for c in sub_cols) if isinstance(s, str) and s]
          for row in complex_table.itertuples(index=False)}
    entities = set(it.ligand) | set(it.receptor)
    return {e: cx.get(e, [e]) for e in entities}


def run_grid(pairs, subunits, datasets, libraries, thresholds=THRESHOLDS):
    """pairs: CellChat Secreted/Contact pairs with ligand_complex/receptor_complex flags."""
    steps, passing = [], []
    for ds_name, ds in datasets.items():
        frac = gene_fractions(ds['X'], ds['n_spots'])
        det = entity_detection(subunits, ds['X'], ds['genes'], frac)
        frac_of = dict(zip(ds['genes'], frac))
        keymaps = {lib.name: genes_by_key(ds['genes'], lib.upper) for lib in libraries}
        for thr in thresholds:
            for rule in COMPLEX_RULES:
                for detection in (['codetect'] if rule == 'single_gene_only' else ['codetect', 'per_subunit']):
                    p = pairs.copy()
                    p['step_a'] = (~p.ligand_complex & ~p.receptor_complex) if rule == 'single_gene_only' else True
                    p['ligand_det'] = p.ligand.map(det[detection])
                    p['receptor_det'] = p.receptor.map(det[detection])
                    p['step_b'] = p.step_a & (p.ligand_det >= thr) & (p.receptor_det >= thr)
                    b = p[p.step_b]
                    excluded = set().union(*[set(subunits[e]) for e in set(b.ligand) | set(b.receptor)]) \
                        if len(b) else set()
                    for lib in libraries:
                        info = {}
                        for rec in sorted(set(b.receptor)):
                            subs = subunits[rec]
                            linked, tfs = lib.receptor(subs)
                            link_ok = (len(linked) == len(subs)) if rule == 'complex_all_subunits_tf' else bool(linked)
                            kept = [g for g in lib.target_symbols(tfs, keymaps[lib.name]) - excluded
                                    if frac_of.get(g, 0) >= thr]
                            info[rec] = (link_ok, len(kept))
                        q = p.copy()
                        q['step_c1'] = q.step_b & q.receptor.map(lambda r: info.get(r, (False, 0))[0])
                        q['n_targets'] = q.receptor.map(lambda r: info.get(r, (False, np.nan))[1])
                        q['step_c'] = q.step_c1 & (q.n_targets >= MIN_TARGETS)
                        meta = dict(dataset=ds_name, threshold=thr, library=lib.name, complex_rule=rule,
                                    detection=detection)
                        for cls in CLASSES:
                            k = q[q.annotation == cls]
                            for step, mask in [('0', np.ones(len(k), bool)), ('a', k.step_a), ('b', k.step_b),
                                               ('c1', k.step_c1), ('c', k.step_c)]:
                                steps.append(dict(**meta, cls=cls, step=step, **sharing(k[np.asarray(mask, bool)])))
                        keep = q[q.step_b]
                        passing.append(keep[['interaction_name', 'annotation', 'pathway_name', 'ligand', 'receptor',
                                             'ligand_complex', 'receptor_complex', 'ligand_det', 'receptor_det',
                                             'step_c1', 'n_targets', 'step_c']].assign(**meta))
    return pd.DataFrame(steps), pd.concat(passing, ignore_index=True)


def markdown_grid(steps, cls, detection='codetect'):
    s = steps[(steps.cls == cls) & (steps.detection == detection)]
    order = {r: i for i, r in enumerate(COMPLEX_RULES)}
    keys = ['dataset', 'threshold', 'library', 'complex_rule']
    lines = ['| 数据集 | 阈值 | 靶基因库 | 复合物规则 | a | b | c1 | c |', '|---|---|---|---|---|---|---|---|']
    groups = sorted(s.groupby(keys), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], order[kv[0][3]]))
    for (ds, thr, lib, rule), g in groups:
        cell = {r.step: f'{r.pairs} / {r.receptors}' for r in g.itertuples()}
        lines.append(f'| {ds} | {thr:.0%} | {lib} | {rule} | {cell["a"]} | {cell["b"]} | {cell["c1"]} | {cell["c"]} |')
    return '\n'.join(lines)


def report(args):
    from cellchat_inventory import load_cellchat
    from evaluate_targets import load_db
    from stage1 import digest, write_json
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    it, cmeta = load_cellchat(args.cellchat)
    complex_table = pd.read_csv(Path(args.cellchat)/'complex.csv', dtype=str)
    pairs = it[it.annotation.isin(CLASSES)].copy()
    subunits = subunit_map(pairs, complex_table)
    sc_rt, sc_tt, sc_meta = load_db(args.scseqcomm)
    cv_rt, cv_tt, cv_meta = load_db(args.cccvelo_db)
    libraries = [Library('scSeqComm', sc_rt, sc_tt, upper=False), Library('CCCvelo', cv_rt, cv_tt, upper=True)]
    vis_h5 = Path(args.visium)/f'{VISIUM_SAMPLE}_filtered_feature_bc_matrix.h5'
    check(digest(vis_h5) == VISIUM_FILES[vis_h5.name], 'Visium matrix changed or missing; run fetch-visium')
    datasets = {'Stereo-seq bin60': load_stereo(args.h5ad, args.matrix), 'Visium posterior': load_visium(vis_h5)}
    steps, passing = run_grid(pairs, subunits, datasets, libraries)
    out.mkdir(parents=True)
    steps.to_csv(out/'steps.csv', index=False)
    passing.to_csv(out/'pairs_passing_b.csv', index=False)
    used = sorted({s for e in set(pairs.ligand) | set(pairs.receptor) for s in subunits[e]})
    dup = {n: [g for g in d['duplicated_symbols'] if g in used] for n, d in datasets.items()}
    write_json(out/'inventory.json', dict(status='complete',
        datasets={n: dict(n_spots=int(d['n_spots']), n_genes=len(d['genes']),
                          duplicated_symbols_summed=len(d['duplicated_symbols']),
                          duplicated_symbols_used_by_cellchat=dup[n]) for n, d in datasets.items()},
        h5ad_sha256=digest(args.h5ad), visium_matrix_sha256=VISIUM_FILES[vis_h5.name],
        cellchat=dict((k, cmeta[k]) for k in ['commit', 'raw_sha256']),
        scseqcomm=sc_meta['source'], cccvelo_db=dict(source=cv_meta['source'], raw_sha256=cv_meta['raw_sha256']),
        thresholds=list(THRESHOLDS), top_tf=TOP_TF, min_targets=MIN_TARGETS, complex_rules=list(COMPLEX_RULES),
        note='Counts only; no panel selected; no OT.'))
    tables = []
    for cls in CLASSES:
        tables.append(f'### {cls}（main detection: co-detection）\n\n' + markdown_grid(steps, cls))
    (out/'grid_tables.md').write_text('\n\n'.join(tables) + '\n')
    print((out/'grid_tables.md').read_text())


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fetch-cccvelo-db')
    f.add_argument('--output', default='data/cccvelo_db'); f.add_argument('--lr-dir', default='data/cccvelo_lr')
    v = sub.add_parser('fetch-visium'); v.add_argument('--output', default='data/visium_mouse_brain_sagittal_posterior')
    r = sub.add_parser('report')
    r.add_argument('--h5ad', default='data/mousebrain_bin60.h5ad'); r.add_argument('--matrix', default='layer:count')
    r.add_argument('--visium', default='data/visium_mouse_brain_sagittal_posterior')
    r.add_argument('--cellchat', default='data/cellchatdb'); r.add_argument('--scseqcomm', default='data/scseqcomm')
    r.add_argument('--cccvelo-db', default='data/cccvelo_db')
    r.add_argument('--output', required=True)
    args = p.parse_args()
    if args.command == 'fetch-cccvelo-db':
        fetch_cccvelo_db(args.output, args.lr_dir)
    elif args.command == 'fetch-visium':
        fetch_visium(args.output)
    else:
        report(args)


if __name__ == '__main__':
    main()
