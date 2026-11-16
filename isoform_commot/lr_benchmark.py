"""PLM-interact on CellChatDB ligand-receptor pairs: data set and negatives (commot-research env).

Plan: reports/plm_lr_benchmark_plan.md (commit 88d1fe1). Evaluation only, never part of OT.

  fetch    pinned downloads into data/lr_benchmark/ (CellChatDB human, CellPhoneDB v5,
           PLM-interact training sets); SHA256 recorded at first download and checked
           afterwards; CellChatDB .rda tables converted to CSV (the mouse tables come
           from the already pinned data/cellchatdb/raw/CellChatDB.mouse.rda)
  genes    gene symbols that need a UniProt entry (single genes and complex subunits
           of the protein interactions), one CSV per species, for
           `proteins.py map --genes ... --species ... --allow-unresolved`
  build    positives, exclusion set E, random (1:10) and degree-preserving (1:1,
           three seeds) negatives for H-M, H-A, H-X, M-M, M-A (plan sections 1-2)

Downloads are data only: CSV via pandas, .rda via the pure-Python `rdata` parser.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd

from research import check

ROOT = Path(__file__).resolve().parent
CELLCHAT_COMMIT = '75253cd0c9e68410e6e721a6d3a0419a1d7e358f'
CELLCHAT_RAW = f'https://raw.githubusercontent.com/jinworks/CellChat/{CELLCHAT_COMMIT}/data/'
MOUSE_RDA = ROOT/'data'/'cellchatdb'/'raw'/'CellChatDB.mouse.rda'
MOUSE_RDA_SHA256 = 'dd59b19646e546b3ed41f394ad2cf64cf44f3d3fcfbb58ac176f0fe1322f9413'
CPDB_TAG, CPDB_COMMIT = 'v5.0.0', '71ffa8a163389913907d40993bec2821a111f64d'
CPDB_RAW = f'https://raw.githubusercontent.com/ventolab/cellphonedb-data/{CPDB_COMMIT}/data/'
CPDB_FILES = ['interaction_input.csv', 'complex_input.csv', 'gene_input.csv', 'protein_input.csv']
HF = 'https://huggingface.co/datasets/danliu1226/{repo}/resolve/{rev}/{path}'
TRAINING = {  # model -> (HF dataset, revision, files); plan 4.7
    'humanV11': ('cross_species_benchmarking', '19dbed25184466fda5341b5ba70c2560ade1bf85',
                 ['train/human.ppi.qrels.seq.train.csv', 'val/human.ppi.qrels.seq.test.csv']),
    'humanV12': ('STRING_V12_TrainingSet', '819f938b1c21b8ec4990bb58c1506cff2537e3c2',
                 ['protein.pairs_human_V12_train.tsv']),
    'LeakageFree': ('Bernett_benchmarking', '5d2ad03baa165c27df32a2eadf066462a2a83073',
                    ['pairs_uniprot_seqs_train.csv', 'pairs_uniprot_seqs_val.csv']),
}
CLASSES = {'Secreted Signaling': 'Secreted', 'Cell-Cell Contact': 'Contact', 'ECM-Receptor': 'ECM'}
NON_PROTEIN = 'Non-protein Signaling'
COFACTOR_COLUMNS = ['agonist', 'antagonist', 'co_A_receptor', 'co_I_receptor']
MAX_AA = 800
RANDOM_RATIO, RANDOM_SEED, DEGREE_SEEDS = 10, 1, (1, 2, 3)
MIN_SWAPS, MAX_ATTEMPTS = 10, 1000          # x number of edges (plan 2.4)
OFFICIAL_MAX_TOKENS = 1603


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def git_state():
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    dirty = subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=ROOT,
                           capture_output=True, text=True).stdout.strip()
    return dict(git_head=head, tracked_changes=bool(dirty))


# ---------------------------------------------------------------- fetch
def fetch_file(url, target, manifest):
    """Download once; afterwards the file must match the SHA256 recorded at first download."""
    from stage1 import digest, download
    key = str(target.relative_to(ROOT))
    known = manifest.get(key, {}).get('sha256')
    download(url, target, expected=known)
    sha = digest(target)
    check(known is None or sha == known, f'{key} changed')
    manifest[key] = dict(url=url, sha256=sha, bytes=target.stat().st_size,
                         retrieved_utc=manifest.get(key, {}).get('retrieved_utc', now()),
                         sha256_source='calculated locally at first download, not author-published')


def rda_tables(path, folder):
    import rdata
    obj = next(iter(rdata.read_rda(path).values()))
    folder.mkdir(parents=True, exist_ok=True)
    interaction = obj['interaction'].copy()
    interaction.columns = [str(c) for c in interaction.columns]
    interaction.to_csv(folder/'interaction.csv', index=False)
    for name in ['complex', 'cofactor']:
        table = obj[name].copy()
        table.index = table.index.astype(str)
        table.columns = [str(c) for c in table.columns]
        table.to_csv(folder/f'{name}.csv', index_label=name)
    obj['geneInfo'].astype(str).to_csv(folder/'gene_info.csv', index=False)
    return dict(n_interactions=len(interaction), annotation_counts=interaction.annotation.value_counts().to_dict())


def fetch(args):
    from stage1 import digest, write_json
    out = ROOT/'data'/'lr_benchmark'
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = out/'sources.json'
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    if 'cellchat' in args.what:
        raw = out/'cellchatdb_human'/'raw'/'CellChatDB.human.rda'
        fetch_file(CELLCHAT_RAW + 'CellChatDB.human.rda', raw, manifest)
        manifest['cellchatdb_human_tables'] = dict(commit=CELLCHAT_COMMIT, **rda_tables(raw, out/'cellchatdb_human'))
        check(digest(MOUSE_RDA) == MOUSE_RDA_SHA256, 'Pinned mouse CellChatDB changed')
        manifest['cellchatdb_mouse_tables'] = dict(commit=CELLCHAT_COMMIT, source=str(MOUSE_RDA.relative_to(ROOT)),
            sha256=MOUSE_RDA_SHA256, **rda_tables(MOUSE_RDA, out/'cellchatdb_mouse'))
    if 'cellphonedb' in args.what:
        for name in CPDB_FILES:
            fetch_file(CPDB_RAW + name, out/'cellphonedb_v5'/name, manifest)
        manifest['cellphonedb_v5'] = dict(tag=CPDB_TAG, commit=CPDB_COMMIT,
                                          license='no license file in the data repository; used only to exclude pairs')
    if 'training' in args.what:
        for model, (repo, rev, files) in TRAINING.items():
            for name in files:
                fetch_file(HF.format(repo=repo, rev=rev, path=name), out/'plm_training'/repo/name, manifest)
            manifest[f'training_{model}'] = dict(dataset=f'danliu1226/{repo}', revision=rev, files=files)
    write_json(manifest_path, manifest)
    print(json.dumps({k: v.get('sha256', v) if isinstance(v, dict) else v for k, v in manifest.items()}, indent=1)[:3000])


# ---------------------------------------------------------------- CellChatDB tables
def clean(values):
    return [str(v).strip() for v in values if pd.notna(v) and str(v).strip() not in ('', 'nan', 'NA')]


def load_cellchat(species):
    folder = ROOT/'data'/'lr_benchmark'/f'cellchatdb_{species}'
    check((folder/'interaction.csv').exists(), f'Run fetch --what cellchat first ({folder})')
    it = pd.read_csv(folder/'interaction.csv', dtype=str, keep_default_na=False)
    cx = pd.read_csv(folder/'complex.csv', dtype=str, keep_default_na=False).set_index('complex')
    cf = pd.read_csv(folder/'cofactor.csv', dtype=str, keep_default_na=False).set_index('cofactor')
    complexes = {name: clean(row) for name, row in cx.iterrows()}
    cofactors = {name: clean(row) for name, row in cf.iterrows()}
    return it, complexes, cofactors


def units(entity, complexes):
    return complexes.get(entity) or [entity]


def protein_genes(species):
    it, complexes, _ = load_cellchat(species)
    prot = it[it.annotation != NON_PROTEIN]
    genes = set()
    for side in ['ligand', 'receptor']:
        for entity in prot[side]:
            genes.update(units(entity, complexes))
    return sorted(genes)


def genes(args):
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    for species in ['human', 'mouse']:
        path = out/f'genes_{species}.csv'
        check(not path.exists(), f'{path} exists')
        pd.DataFrame(dict(symbol=protein_genes(species))).to_csv(path, index=False)
        print(f'{species}: {len(pd.read_csv(path))} genes -> {path}')


# ---------------------------------------------------------------- mapping
def load_map(folder):
    folder = Path(folder)
    meta = json.loads((folder/'proteins.json').read_text())
    check(meta.get('status') == 'complete', f'{folder}: incomplete protein map')
    table = pd.read_csv(folder/'protein_map.csv', dtype=str)
    table['length'] = table.length.astype(int)
    return table.set_index('lr_symbol'), meta


# ---------------------------------------------------------------- exclusion set
class Exclusion:
    """Unordered pairs that may not be negatives: symbol level (upper case) and accession level."""

    def __init__(self):
        self.sym, self.acc = {}, {}

    @staticmethod
    def _add(store, a, b, source):
        if a and b:
            store.setdefault(frozenset((a, b)), set()).add(source)

    def add_symbols(self, left, right, source):
        for a in left:
            for b in right:
                self._add(self.sym, a.upper(), b.upper(), source)

    def add_accessions(self, left, right, source):
        for a in left:
            for b in right:
                self._add(self.acc, a, b, source)

    def sources(self, sym_a, sym_b, acc_a=None, acc_b=None):
        found = set(self.sym.get(frozenset((sym_a.upper(), sym_b.upper())), ()))
        if acc_a and acc_b:
            found |= self.acc.get(frozenset((acc_a, acc_b)), set())
        return found


def cellchat_exclusions(excl, species):
    it, complexes, cofactors = load_cellchat(species)
    missing_cofactors = set()
    for r in it.to_dict('records'):
        left = set(units(r['ligand'], complexes)) | set(clean(r['ligand.symbol'].split(',')))
        right = set(units(r['receptor'], complexes)) | set(clean(str(r['receptor.symbol']).split(',')))
        excl.add_symbols(left, right, f'cellchat_{species}')
        for col in COFACTOR_COLUMNS:
            name = r[col].strip()
            if name:
                if name not in cofactors:
                    missing_cofactors.add(name)
                excl.add_symbols(left, cofactors.get(name, []), f'cellchat_cofactor_{species}')
    return dict(rows=len(it), missing_cofactor_names=sorted(missing_cofactors))


def cellphonedb_exclusions(excl):
    folder = ROOT/'data'/'lr_benchmark'/'cellphonedb_v5'
    every = pd.read_csv(folder/'interaction_input.csv', dtype=str, keep_default_na=False)
    inter = every[every.is_ppi == 'True']  # plan 2.1: protein-protein only; small-molecule partners ignored
    cx = pd.read_csv(folder/'complex_input.csv', dtype=str, keep_default_na=False)
    gene = pd.read_csv(folder/'gene_input.csv', dtype=str, keep_default_na=False)
    protein = pd.read_csv(folder/'protein_input.csv', dtype=str, keep_default_na=False)
    sub_cols = [c for c in cx.columns if c.startswith('uniprot_')]
    complexes = {r['complex_name']: clean(r[c] for c in sub_cols) for r in cx.to_dict('records')}
    accessions = set(protein.uniprot)
    symbols = {}
    for r in gene.to_dict('records'):
        for col in ['gene_name', 'hgnc_symbol']:
            if r.get(col):
                symbols.setdefault(r['uniprot'], set()).add(r[col])
    skipped, used = Counter(), 0
    for r in inter.to_dict('records'):
        sides = []
        for partner in [r['partner_a'], r['partner_b']]:
            if partner in complexes:
                sides.append(complexes[partner])
            elif partner in accessions:
                sides.append([partner])
            else:
                sides.append(None)
                skipped[partner] += 1
        if sides[0] is None or sides[1] is None:
            continue
        used += 1
        excl.add_accessions(sides[0], sides[1], 'cellphonedb_v5')
        sym = [set().union(*[symbols.get(a, set()) for a in side]) for side in sides]
        excl.add_symbols(sym[0], sym[1], 'cellphonedb_v5')
    return dict(interactions=len(every), is_ppi_true=len(inter), used_protein_interactions=used,
                skipped_non_protein_partners=int(sum(skipped.values())), n_skipped_partner_names=len(skipped),
                columns=inter.columns.tolist()[:12])


# ---------------------------------------------------------------- positives
def positives(species, pmap, unresolved):
    """Unordered entity pairs of protein interactions, merged; first row gives the orientation."""
    it, complexes, _ = load_cellchat(species)
    prot = it[it.annotation != NON_PROTEIN]
    merged = {}
    symbol_mismatch = []
    for r in prot.to_dict('records'):
        lig, rec = r['ligand'], r['receptor']
        single = lig not in complexes and rec not in complexes
        if single and (r['ligand.symbol'].strip() != lig or str(r['receptor.symbol']).strip() != rec):
            symbol_mismatch.append(r['interaction_name'])
        key = frozenset((lig, rec))
        if key not in merged:
            merged[key] = dict(ligand=lig, receptor=rec, single=single, classes=set(), pathways=set(),
                               evidence=set(), interactions=[], both_orientations=False)
        m = merged[key]
        m['classes'].add(CLASSES[r['annotation']])
        m['pathways'].add(r['pathway_name'])
        m['evidence'].add(r['evidence'])
        m['interactions'].append(r['interaction_name'])
        if (lig, rec) != (m['ligand'], m['receptor']):
            m['both_orientations'] = True
    rows = []
    for m in merged.values():
        lu, ru = units(m['ligand'], complexes), units(m['receptor'], complexes)
        missing = sorted({g for g in lu + ru if g in unresolved or g not in pmap.index})
        row = dict(ligand=m['ligand'], receptor=m['receptor'], single=m['single'],
                   homophilic=m['ligand'] == m['receptor'], classes=';'.join(sorted(m['classes'])),
                   pathways=';'.join(sorted(m['pathways'])), evidence=' | '.join(sorted(m['evidence'])),
                   interactions=';'.join(m['interactions']), both_orientations=m['both_orientations'],
                   ligand_subunits=';'.join(lu), receptor_subunits=';'.join(ru),
                   unmapped=';'.join(missing))
        if not missing:
            la = [pmap.loc[g, 'accession'] for g in lu]
            ra = [pmap.loc[g, 'accession'] for g in ru]
            ll = [int(pmap.loc[g, 'length']) for g in lu]
            rl = [int(pmap.loc[g, 'length']) for g in ru]
            row.update(ligand_accessions=';'.join(la), receptor_accessions=';'.join(ra),
                       ligand_lengths=';'.join(map(str, ll)), receptor_lengths=';'.join(map(str, rl)),
                       max_len=max(ll + rl),
                       n_subunit_pairs=sum(a != b for a in la for b in ra))
        rows.append(row)
    table = pd.DataFrame(rows).sort_values(['ligand', 'receptor']).reset_index(drop=True)
    return table, symbol_mismatch


# ---------------------------------------------------------------- negatives
class Entities:
    """Entity -> subunit symbols / accessions / lengths, and entity-pair exclusion."""

    def __init__(self, species, pmap, complexes, excl):
        self.species, self.pmap, self.complexes, self.excl = species, pmap, complexes, excl
        self.cache = {}

    def subunits(self, entity):
        return units(entity, self.complexes)

    def accessions(self, entity):
        return [self.pmap.loc[g, 'accession'] for g in self.subunits(entity)]

    def lengths(self, entity):
        return [int(self.pmap.loc[g, 'length']) for g in self.subunits(entity)]

    def sources(self, a, b):
        key = frozenset((a, b))
        if key not in self.cache:
            found = set()
            human = self.species == 'human'
            for ga in self.subunits(a):
                for gb in self.subunits(b):
                    if ga.upper() == gb.upper():
                        continue
                    acc_a = self.pmap.loc[ga, 'accession'] if human else None
                    acc_b = self.pmap.loc[gb, 'accession'] if human else None
                    found |= self.excl.sources(ga, gb, acc_a, acc_b)
            self.cache[key] = found
        return self.cache[key]

    def excluded(self, a, b):
        return bool(self.sources(a, b))

    def has_distinct_pair(self, a, b):
        return any(x != y for x in self.accessions(a) for y in self.accessions(b))


def candidates(pos, ent):
    """Ligand pool x receptor pool, L != R, not in E, unordered unique (plan 2.2)."""
    ligands, receptors = sorted(set(pos.ligand)), sorted(set(pos.receptor))
    positive_keys = {frozenset(k) for k in zip(pos.ligand, pos.receptor)}
    seen, out, removed, annotated = set(), [], Counter(), []
    for l in ligands:
        for r in receptors:
            key = frozenset((l, r))
            if l == r or key in seen:
                continue
            seen.add(key)
            src = ent.sources(l, r)
            if src:
                for s in src:
                    removed[s] += 1
                annotated.append(dict(ligand=l, receptor=r, positive=key in positive_keys,
                                      sources=';'.join(sorted(src))))
                continue
            if not ent.has_distinct_pair(l, r):
                removed['same_protein'] += 1
                continue
            check(key not in positive_keys, f'positive {l}-{r} not in E')
            out.append((l, r))
    return out, ligands, receptors, removed, pd.DataFrame(annotated)


def degree_negatives(pos, ent, seed):
    """Double-edge swaps on the ligand-receptor bipartite graph (plan 2.4)."""
    rng = np.random.default_rng(seed)
    edges = list(zip(pos.ligand, pos.receptor))
    m = len(edges)
    present = {frozenset(e) for e in edges}
    bad = [ent.excluded(l, r) for l, r in edges]
    n_bad, successes, attempts = sum(bad), 0, 0
    while (n_bad > 0 or successes < MIN_SWAPS*m) and attempts < MAX_ATTEMPTS*m:
        attempts += 1
        i, j = (int(x) for x in rng.integers(m, size=2))
        if i == j:
            continue
        (l1, r1), (l2, r2) = edges[i], edges[j]
        a, b = (l1, r2), (l2, r1)
        if l1 == r2 or l2 == r1:
            continue
        ka, kb = frozenset(a), frozenset(b)
        if ka == kb or ka in present or kb in present:
            continue
        if ent.excluded(*a) or ent.excluded(*b) or not ent.has_distinct_pair(*a) or not ent.has_distinct_pair(*b):
            continue
        present -= {frozenset(edges[i]), frozenset(edges[j])}
        present |= {ka, kb}
        n_bad -= bad[i] + bad[j]
        edges[i], edges[j], bad[i], bad[j] = a, b, False, False
        successes += 1
    kept = [e for e, x in zip(edges, bad) if not x]
    dropped = [e for e, x in zip(edges, bad) if x]
    deg_pos_l, deg_pos_r = Counter(pos.ligand), Counter(pos.receptor)
    deg_l, deg_r = Counter(l for l, _ in kept), Counter(r for _, r in kept)
    deviation = {f'L:{g}': deg_l[g] - n for g, n in deg_pos_l.items() if deg_l[g] != n}
    deviation.update({f'R:{g}': deg_r[g] - n for g, n in deg_pos_r.items() if deg_r[g] != n})
    stats = dict(seed=seed, edges=m, successful_swaps=successes, attempts=attempts,
                 swaps_per_edge=round(successes/m, 2), dropped_still_excluded=len(dropped),
                 degree_deviation=deviation,
                 unchanged_from_start=0)
    return kept, stats


def describe(pairs, ent, pos_ligand_classes, kind):
    rows = []
    for l, r in pairs:
        la, ra, ll, rl = ent.accessions(l), ent.accessions(r), ent.lengths(l), ent.lengths(r)
        rows.append(dict(ligand=l, receptor=r, ligand_accessions=';'.join(la), receptor_accessions=';'.join(ra),
                         ligand_lengths=';'.join(map(str, ll)), receptor_lengths=';'.join(map(str, rl)),
                         max_len=max(ll + rl), both_le_800=max(ll + rl) <= MAX_AA,
                         n_subunit_pairs=sum(a != b for a in la for b in ra),
                         ligand_classes=pos_ligand_classes.get(l, ''), negative=kind))
    return pd.DataFrame(rows)


def ligand_classes(pos):
    out = {}
    for l, classes in zip(pos.ligand, pos.classes):
        out.setdefault(l, set()).update(classes.split(';'))
    return {l: ';'.join(sorted(c)) for l, c in out.items()}


def build_set(name, pos, ent, out):
    folder = out/name
    folder.mkdir(parents=True)
    pos = pos.sort_values(['ligand', 'receptor']).reset_index(drop=True)
    pos.to_csv(folder/'positives.csv', index=False)
    classes = ligand_classes(pos)
    cands, ligands, receptors, removed, annotated = candidates(pos, ent)
    annotated.to_csv(folder/'pool_annotated.csv', index=False)  # E inside the pool (ranking rule, plan 4.3)
    rng = np.random.default_rng(RANDOM_SEED)
    n_random = RANDOM_RATIO * len(pos)
    take = rng.choice(len(cands), size=min(n_random, len(cands)), replace=False)
    random_pairs = [cands[i] for i in take]
    describe(random_pairs, ent, classes, 'random').to_csv(folder/f'negatives_random_s{RANDOM_SEED}.csv', index=False)
    degree = {}
    for seed in DEGREE_SEEDS:
        kept, stats = degree_negatives(pos, ent, seed)
        start = {frozenset(e) for e in zip(pos.ligand, pos.receptor)}
        stats['unchanged_from_start'] = sum(frozenset(e) in start for e in kept)
        check(stats['unchanged_from_start'] == 0, 'a positive survived as a negative')
        describe(kept, ent, classes, f'degree_s{seed}').to_csv(folder/f'negatives_degree_s{seed}.csv', index=False)
        degree[f's{seed}'] = stats
    grid = sum(1 for l in ligands for r in receptors if l != r)
    summary = dict(set=name, positives=len(pos), ligand_pool=len(ligands), receptor_pool=len(receptors),
                   pool_pairs_excluding_self=grid, candidates_after_exclusion=len(cands),
                   removed_by_source=dict(removed), random=dict(seed=RANDOM_SEED, ratio=RANDOM_RATIO,
                   requested=n_random, drawn=len(random_pairs)), degree=degree)
    return summary


def class_counts(pos):
    counts = Counter()
    for c in pos.classes:
        for k in c.split(';'):
            counts[k] += 1
    return {k: counts.get(k, 0) for k in CLASSES.values()} | dict(multi_class=int((pos.classes.str.count(';') > 0).sum()))


def build(args):
    from stage1 import digest, write_json
    out = Path(args.output)
    check(out.is_dir() and not (out/'build.json').exists(), f'{out} must exist (genes step) and have no build yet')
    excl = Exclusion()
    sources = dict(cellchat_human=cellchat_exclusions(excl, 'human'),
                   cellchat_mouse=cellchat_exclusions(excl, 'mouse'),
                   cellphonedb_v5=cellphonedb_exclusions(excl))
    counts, sets, mismatches = [], {}, {}
    for species, prefix, map_dir in [('human', 'H', args.human_map), ('mouse', 'M', args.mouse_map)]:
        pmap, meta = load_map(map_dir)
        unresolved = set(meta.get('unresolved') or [])
        pos, mismatch = positives(species, pmap, unresolved)
        mismatches[species] = mismatch
        pos.to_csv(out/f'positives_all_{species}.csv', index=False)
        _, complexes, _ = load_cellchat(species)
        ent = Entities(species, pmap, complexes, excl)
        mapped = pos[pos.unmapped == '']
        single = mapped[mapped.single & ~mapped.homophilic]
        definitions = {
            f'{prefix}-A': single,
            f'{prefix}-M': single[single.max_len <= MAX_AA],
        }
        if species == 'human':
            x = mapped[~mapped.homophilic]
            x = x[[ent.has_distinct_pair(l, r) for l, r in zip(x.ligand, x.receptor)]]
            definitions['H-X'] = x
        mapped[mapped.homophilic & mapped.single].to_csv(out/f'homophilic_{species}.csv', index=False)
        for name, table in definitions.items():
            sets[name] = build_set(name, table, ent, out)
            accessions = set()
            for col in ['ligand_accessions', 'receptor_accessions']:
                for v in table[col]:
                    accessions.update(v.split(';'))
            counts.append(dict(set=name, species=species, positives=len(table), **class_counts(table),
                               complex_pairs=int((~table.single).sum()),
                               proteins=len(accessions),
                               both_le_800=int((table.max_len <= MAX_AA).sum()),
                               frac_both_le_800=round(float((table.max_len <= MAX_AA).mean()), 4),
                               random_negatives=sets[name]['random']['drawn'],
                               **{f'degree_negatives_s{s}': len(pd.read_csv(out/name/f'negatives_degree_s{s}.csv'))
                                  for s in DEGREE_SEEDS},
                               ranking_grid=sets[name]['pool_pairs_excluding_self']
                               if name.endswith('-M') else None))
        all_protein = len(pos)
        counts.append(dict(set=f'{prefix} (all protein pairs, before rules)', species=species,
                           positives=all_protein, **class_counts(pos),
                           complex_pairs=int((~pos.single).sum()),
                           homophilic_single=int((pos.single & pos.homophilic).sum()),
                           dropped_unmapped=int((pos.unmapped != '').sum()),
                           dropped_unmapped_single=int(((pos.unmapped != '') & pos.single).sum()),
                           unresolved_genes=len(unresolved)))
    pd.DataFrame(counts).to_csv(out/'counts.csv', index=False)
    import shutil  # provenance copies: protein maps (no sequences), unresolved genes, download hashes
    (out/'sources').mkdir()
    for species, map_dir in [('human', Path(args.human_map)), ('mouse', Path(args.mouse_map))]:
        shutil.copy(map_dir/'protein_map.csv', out/'sources'/f'protein_map_{species}.csv')
        shutil.copy(map_dir/'proteins.json', out/'sources'/f'proteins_{species}.json')
        candidates_file = Path(str(map_dir) + '_candidates.csv')
        if candidates_file.exists():
            shutil.copy(candidates_file, out/'sources'/f'unresolved_{species}.csv')
    shutil.copy(ROOT/'data'/'lr_benchmark'/'sources.json', out/'sources'/'downloads.json')
    files = sorted(p for p in out.rglob('*.csv'))
    write_json(out/'build.json', dict(status='complete', created_utc=now(), **git_state(),
        plan='reports/plm_lr_benchmark_plan.md', human_map=str(args.human_map), mouse_map=str(args.mouse_map),
        human_map_fasta_sha256=load_map(args.human_map)[1]['fasta_sha256'],
        mouse_map_fasta_sha256=load_map(args.mouse_map)[1]['fasta_sha256'],
        exclusion=dict(symbol_pairs=len(excl.sym), accession_pairs=len(excl.acc), sources=sources),
        single_gene_symbol_column_mismatches=mismatches, sets=sets,
        rules=dict(random_ratio=RANDOM_RATIO, random_seed=RANDOM_SEED, degree_seeds=list(DEGREE_SEEDS),
                   min_swaps_per_edge=MIN_SWAPS, max_attempts_per_edge=MAX_ATTEMPTS, max_aa=MAX_AA),
        files={str(p.relative_to(out)): digest(p) for p in files}))
    print(pd.DataFrame(counts).to_string(index=False))


# ---------------------------------------------------------------- scoring jobs
def set_pairs(folder, grid=False):
    """Unordered accession pairs of one analysis set (positives, all negatives, optional ranking grid)."""
    folder = Path(folder)
    pairs = set()
    tables = [pd.read_csv(folder/'positives.csv', dtype=str)]
    tables += [pd.read_csv(p, dtype=str) for p in sorted(folder.glob('negatives_*.csv'))]
    for t in tables:
        for la, ra in zip(t.ligand_accessions, t.receptor_accessions):
            pairs.update(frozenset((a, b)) for a in la.split(';') for b in ra.split(';') if a != b)
    if grid:
        pos = tables[0]
        ligands = dict(zip(pos.ligand, pos.ligand_accessions))
        receptors = dict(zip(pos.receptor, pos.receptor_accessions))
        pairs.update(frozenset((ligands[l], receptors[r])) for l in ligands for r in receptors
                     if ligands[l] != receptors[r])
    return pairs


def ecd_pairs(folder):
    """H-M positives and negatives only (plan 3.5: no ranking on extracellular regions)."""
    return set_pairs(folder, grid=False)


def write_jobs(name, pairs, sequences, out, chunk):
    folder = out/name
    folder.mkdir(parents=True)
    rows = []
    for key in sorted(tuple(sorted(k)) for k in pairs):
        a, b = key
        for q, t in [(a, b), (b, a)]:
            rows.append(dict(pair_id=f'{q}__{t}', query_accession=q, text_accession=t,
                             query=sequences[q], text=sequences[t],
                             n_tokens=len(sequences[q]) + len(sequences[t]) + 3))
    jobs = pd.DataFrame(rows).sort_values(['n_tokens', 'pair_id']).reset_index(drop=True)
    jobs['chunk'] = np.arange(len(jobs)) // chunk
    for c, part in jobs.groupby('chunk'):
        part[['pair_id', 'query', 'text']].to_csv(folder/f'chunk_{c:03d}.csv', index=False)
    jobs.drop(columns=['query', 'text']).to_csv(folder/'index.csv', index=False)
    n = jobs.n_tokens
    return dict(unordered_pairs=len(pairs), jobs=len(jobs), chunks=int(jobs.chunk.max()) + 1,
                over_1603_tokens=int((n > OFFICIAL_MAX_TOKENS).sum()), tokens_median=int(n.median()),
                tokens_max=int(n.max()), tokens_sum=int(n.sum()), tokens_sq_sum=float((n.astype(float)**2).sum()))


def jobs(args):
    from proteins import region_sequence
    from stage1 import write_json
    data, out = Path(args.dataset), Path(args.output)
    check(not out.exists(), f'{out} exists')
    out.mkdir(parents=True)
    summary = {}
    for species, map_dir, sets in [('human', args.human_map, [('H-M', True), ('H-A', False), ('H-X', False)]),
                                    ('mouse', args.mouse_map, [('M-M', True), ('M-A', False)])]:
        entries = {p.stem: json.loads(p.read_text()) for p in (Path(map_dir)/'uniprot').glob('*.json')}
        full = {a: e['sequence']['value'] for a, e in entries.items()}
        pairs = set()
        for name, grid in sets:
            pairs |= set_pairs(data/name, grid)
        summary[f'{species}_full'] = write_jobs(f'{species}_full', pairs, full, out, args.chunk)
        if species == 'human':
            ecd, refused = {}, {}
            for a, e in entries.items():
                try:
                    ecd[a] = region_sequence(e, 'extracellular')[0]
                except ValueError as exc:
                    refused[a] = str(exc)
            hm = ecd_pairs(data/'H-M')
            usable = {k for k in hm if all(a in ecd for a in k)}
            summary['human_ecd'] = write_jobs('human_ecd', usable, ecd, out, args.chunk)
            summary['human_ecd'].update(hm_pairs=len(hm), hm_pairs_without_ecd=len(hm) - len(usable),
                                        proteins_refused=len(refused))
            pd.Series(refused, name='reason').rename_axis('accession').to_csv(out/'ecd_refused.csv')
    write_json(out/'jobs.json', dict(status='complete', created_utc=now(), **git_state(), dataset=str(data),
                                     chunk=args.chunk, summary=summary))
    print(json.dumps(summary, indent=1))


# ---------------------------------------------------------------- scoring workers
CHECKPOINTS = {  # model -> (local file, HF name, revision); plan 3.1
    'humanV11': ('models/plm_650m_humanV11/pytorch_model.bin', 'danliu1226/PLM-interact-650M-humanV11',
                 'e86e392dec13dd0c23252c94947b04a7a9821b0e'),
    'humanV12': ('models/plm_650m_humanV12/pytorch_model.bin', 'danliu1226/PLM-interact-650M-humanV12',
                 '2d0e7fb3deef94efa3f217daccc79092a6a42ceb'),
    'LeakageFree': ('models/plm_650m_leakagefree/pytorch_model.bin',
                    'danliu1226/PLM-interact-650M-Leakage-Free-Dataset', 'd995319a6e53510bfbb4eff6ab0292fedf2ec4db'),
}
SCORING = [  # (jobset, model, scorers); plan 3.5
    ('human_full', 'humanV11', ['ours', 'official']), ('human_full', 'humanV12', ['ours', 'official']),
    ('human_full', 'LeakageFree', ['ours', 'official']), ('mouse_full', 'humanV11', ['ours', 'official']),
    ('human_ecd', 'humanV11', ['ours']),
]
PLM_PYTHON = '/home/liu/.conda/envs/plminteract/bin/python'


def tasks(runs):
    out = []
    for jobset, model, scorers in SCORING:
        for chunk in sorted((runs/'jobs'/jobset).glob('chunk_*.csv')):
            for scorer in scorers:
                out.append((jobset, model, scorer, chunk))
    return out


def work(args):
    """Worker w of n: every n-th task; finished tasks (score file + meta) are skipped, so it can resume."""
    runs = Path(args.runs)
    ckpt = {m: (str(ROOT/f), n, r) for m, (f, n, r) in CHECKPOINTS.items()}
    for i, (jobset, model, scorer, chunk) in enumerate(tasks(runs)):
        if i % args.of != args.worker:
            continue
        target = runs/'scores'/jobset/model/scorer/chunk.name
        if Path(str(target) + '.meta.json').exists():
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        log = target.with_suffix('.log')
        path, name, rev = ckpt[model]
        common = ['--checkpoint', path, '--checkpoint-name', name, '--checkpoint-revision', rev]
        print(f'{now()} worker {args.worker}: {jobset} {model} {scorer} {chunk.name}', flush=True)
        with log.open('w') as stream:
            if scorer == 'ours':
                if target.exists():  # interrupted run: plm_score.py refuses to overwrite
                    target.rename(target.with_suffix(f'.incomplete_{datetime.now():%H%M%S}'))
                cmd = [PLM_PYTHON, 'plm_score.py', 'score', '--jobs', str(chunk), '--base-model',
                       'models/esm2_t33_650M_UR50D', '--allow-long', '--output', str(target)] + common
                subprocess.run(cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, check=True)
            else:
                work_dir = runs/'official'/jobset/model/chunk.stem
                if work_dir.exists():  # interrupted run: start the chunk again in a fresh directory
                    work_dir.rename(work_dir.with_name(f'{work_dir.name}.incomplete_{datetime.now():%H%M%S}'))
                for cmd in [[PLM_PYTHON, 'plm_official.py', 'prepare', '--jobs', str(chunk), '--dir', str(work_dir)],
                            [PLM_PYTHON, 'plm_official.py', 'run', '--dir', str(work_dir)] + common,
                            [PLM_PYTHON, 'plm_official.py', 'collect', '--dir', str(work_dir), '--output', str(target)]]:
                    subprocess.run(cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, check=True)
    print(f'{now()} worker {args.worker}: done', flush=True)


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fetch')
    f.add_argument('--what', nargs='+', choices=['cellchat', 'cellphonedb', 'training'],
                   default=['cellchat', 'cellphonedb', 'training'])
    g = sub.add_parser('genes')
    g.add_argument('--output', default='reports/plm_lr_benchmark/dataset')
    b = sub.add_parser('build')
    b.add_argument('--human-map', required=True)
    b.add_argument('--mouse-map', required=True)
    b.add_argument('--output', default='reports/plm_lr_benchmark/dataset')
    j = sub.add_parser('jobs')
    j.add_argument('--dataset', default='reports/plm_lr_benchmark/dataset')
    j.add_argument('--human-map', required=True)
    j.add_argument('--mouse-map', required=True)
    j.add_argument('--output', default='plm_lr_benchmark_runs/jobs')
    j.add_argument('--chunk', type=int, default=5000)
    w = sub.add_parser('work')
    w.add_argument('--runs', default='plm_lr_benchmark_runs')
    w.add_argument('--worker', type=int, required=True)
    w.add_argument('--of', type=int, default=3)
    return p


if __name__ == '__main__':
    a = parser().parse_args()
    {'fetch': fetch, 'genes': genes, 'build': build, 'jobs': jobs, 'work': work}[a.command](a)
