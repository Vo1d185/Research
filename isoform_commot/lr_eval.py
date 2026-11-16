"""PLM-interact LR benchmark: metrics, strata, ranking, hub diagnostic, consistency, training overlap.

Plan: reports/plm_lr_benchmark_plan.md, sections 4-5 (commot-research env). Inputs: the data set
(reports/plm_lr_benchmark/dataset, lr_benchmark.py build), the scoring jobs (lr_benchmark.py jobs)
and score files laid out as

  <runs>/scores/<jobset>/<model>/ours/chunk_XXX.csv (+ .meta.json, plm_score.py)
  <runs>/scores/<jobset>/<model>/official/chunk_XXX.csv (+ .meta.json, plm_official.py collect)

  collect   chunk files -> <runs>/scores/<jobset>__<model>__<scorer>.csv.gz, checking that every
            job was scored exactly once from the recorded input
  evaluate  everything in plan section 4, and the verdict of section 5 -> <output>/
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

from research import check

MODELS = ['humanV11', 'humanV12', 'LeakageFree']
MAIN_MODEL, MAIN_SET = 'humanV11', 'H-M'
N_BOOT, BOOT_SEED, MAX_AA, TOL = 2000, 1, 800, 1e-3
CRITERIA = dict(random_auroc=0.85, degree_auroc=0.75, top5=0.50)        # plan 5 (pending confirmation)
CLASSES = ['Secreted', 'Contact', 'ECM']
SCOPE = {  # set -> (jobset, models); plan 3.5
    'H-M': ('human_full', MODELS), 'H-A': ('human_full', MODELS), 'H-X': ('human_full', MODELS),
    'M-M': ('mouse_full', ['humanV11']), 'M-A': ('mouse_full', ['humanV11']),
}
TRAINING = {
    'humanV11': ['cross_species_benchmarking/train/human.ppi.qrels.seq.train.csv',
                 'cross_species_benchmarking/val/human.ppi.qrels.seq.test.csv'],
    'humanV12': ['STRING_V12_TrainingSet/protein.pairs_human_V12_train.tsv'],
    'LeakageFree': ['Bernett_benchmarking/pairs_uniprot_seqs_train.csv',
                    'Bernett_benchmarking/pairs_uniprot_seqs_val.csv'],
}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False, default=str))


# ---------------------------------------------------------------- collect
def collect(args):
    runs = Path(args.runs)
    for jobset_dir in sorted((runs/'jobs').iterdir()):
        if not jobset_dir.is_dir():
            continue
        index = pd.read_csv(jobset_dir/'index.csv', dtype={'pair_id': str})
        for scorer_dir in sorted((runs/'scores'/jobset_dir.name).glob('*/*')):
            model, scorer = scorer_dir.parent.name, scorer_dir.name
            out = runs/'scores'/f'{jobset_dir.name}__{model}__{scorer}.csv.gz'
            if out.exists():
                continue
            parts = []
            for chunk in sorted(jobset_dir.glob('chunk_*.csv')):
                path = scorer_dir/chunk.name
                check(path.exists(), f'missing {path}')
                meta = json.loads(Path(str(path) + '.meta.json').read_text())
                check(meta.get('status') == 'complete', f'{path} incomplete')
                source = meta.get('input_sha256') or meta['prepare']['jobs_sha256']
                check(source == digest(chunk), f'{path} was not scored from {chunk}')
                parts.append(pd.read_csv(path, dtype={'pair_id': str}))
            scores = pd.concat(parts, ignore_index=True)
            check(scores.pair_id.is_unique and set(scores.pair_id) == set(index.pair_id),
                  f'{scorer_dir}: scores must cover every job exactly once')
            s = scores.score.to_numpy(float)
            check(np.isfinite(s).all() and ((s >= 0) & (s <= 1)).all(), 'scores outside [0, 1]')
            keep = [c for c in ['pair_id', 'score', 'n_tokens', 'truncated', 'n_tokens_used'] if c in scores]
            scores[keep].to_csv(out, index=False)
            print(f'{out}: {len(scores)} jobs')


# ---------------------------------------------------------------- metrics
def auroc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if not len(pos) or not len(neg):
        return np.nan
    ranks = rankdata(np.concatenate([pos, neg]))           # average ranks: ties count 1/2
    return float((ranks[:len(pos)].sum() - len(pos)*(len(pos) + 1)/2) / (len(pos)*len(neg)))


def aupr(pos, neg):
    """Average precision with tied scores grouped into one threshold (as scikit-learn)."""
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if not len(pos) or not len(neg):
        return np.nan
    scores = np.concatenate([pos, neg])
    labels = np.concatenate([np.ones(len(pos)), np.zeros(len(neg))])
    order = np.argsort(-scores, kind='stable')
    scores, labels = scores[order], labels[order]
    last = np.r_[np.flatnonzero(np.diff(scores)), len(scores) - 1]   # last index of each tie group
    tp = np.cumsum(labels)[last]
    fp = (last + 1) - tp
    precision, recall = tp/(tp + fp), tp/len(pos)
    return float(np.sum(np.diff(np.r_[0, recall]) * precision))


def boot_ci(stat, pos, neg, rng):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if not len(pos) or not len(neg):
        return np.nan, np.nan
    values = [stat(pos[rng.integers(len(pos), size=len(pos))], neg[rng.integers(len(neg), size=len(neg))])
              for _ in range(N_BOOT)]
    return tuple(float(x) for x in np.percentile(values, [2.5, 97.5]))


def metric_row(pos, neg, **keys):
    rng = np.random.default_rng(BOOT_SEED)
    a_lo, a_hi = boot_ci(auroc, pos, neg, rng)
    p_lo, p_hi = boot_ci(aupr, pos, neg, rng)
    return dict(**keys, n_pos=len(pos), n_neg=len(neg), auroc=auroc(pos, neg), auroc_lo=a_lo, auroc_hi=a_hi,
                aupr=aupr(pos, neg), aupr_lo=p_lo, aupr_hi=p_hi,
                aupr_baseline=len(pos)/(len(pos) + len(neg)) if len(pos) + len(neg) else np.nan,
                pos_gt_05=float(np.mean(np.asarray(pos) > 0.5)) if len(pos) else np.nan,
                neg_gt_05=float(np.mean(np.asarray(neg) > 0.5)) if len(neg) else np.nan)


# ---------------------------------------------------------------- scores
class Scores:
    """Ordered probabilities -> unordered pair score (mean of both orders) -> entity pair (max)."""

    def __init__(self, path):
        table = pd.read_csv(path, dtype={'pair_id': str})
        self.p = dict(zip(table.pair_id, table.score.astype(float)))
        self.table = table

    def pair(self, a, b):
        return (self.p[f'{a}__{b}'] + self.p[f'{b}__{a}']) / 2

    def entity(self, ligand_accessions, receptor_accessions):
        values = [self.pair(a, b) for a in ligand_accessions.split(';') for b in receptor_accessions.split(';')
                  if a != b]
        return max(values)


def score_table(folder, scores, extra=None):
    """Positives and every negative file of one set, with scores (extra: other scorers/regions)."""
    folder = Path(folder)
    pos = pd.read_csv(folder/'positives.csv', dtype=str)
    pos['negative'] = 'positive'
    tables = [pos] + [pd.read_csv(p, dtype=str) for p in sorted(folder.glob('negatives_*.csv'))]
    t = pd.concat(tables, ignore_index=True)
    t['label'] = (t.negative == 'positive').astype(int)
    t['max_len'] = t.max_len.astype(int)
    t['complex'] = t.ligand_accessions.str.contains(';') | t.receptor_accessions.str.contains(';')
    t['score'] = [scores.entity(l, r) for l, r in zip(t.ligand_accessions, t.receptor_accessions)]
    for name, other in (extra or {}).items():
        values = []
        for l, r in zip(t.ligand_accessions, t.receptor_accessions):
            try:
                values.append(other.entity(l, r))
            except KeyError:
                values.append(np.nan)
        t[name] = values
    return t


def class_masks(t, k):
    pos = (t.label == 1) & t.classes.fillna('').str.split(';').apply(lambda c: k in c)
    neg = (t.label == 0) & t.ligand_classes.fillna('').str.split(';').apply(lambda c: k in c)
    return pos, neg


def set_metrics(t, keys, column='score'):
    rows = []
    negs = sorted(n for n in t.negative.unique() if n != 'positive')
    for neg_type in negs:
        sub = t[(t.label == 1) | (t.negative == neg_type)]
        p, n = sub[sub.label == 1][column], sub[sub.label == 0][column]
        rows.append(metric_row(p.dropna(), n.dropna(), **keys, negatives=neg_type, stratum='all'))
        for k in CLASSES:
            pm, nm = class_masks(sub, k)
            rows.append(metric_row(sub[pm][column].dropna(), sub[nm][column].dropna(), **keys,
                                   negatives=neg_type, stratum=f'class:{k}'))
        for name, cond in [('both<=800', sub.max_len <= MAX_AA), ('other>800', sub.max_len > MAX_AA)]:
            if cond.all() or not cond.any():
                continue
            s2 = sub[cond]
            rows.append(metric_row(s2[s2.label == 1][column].dropna(), s2[s2.label == 0][column].dropna(), **keys,
                                   negatives=neg_type, stratum=f'length:{name}'))
        if sub.complex.any():
            s2 = sub[sub.complex]
            rows.append(metric_row(s2[s2.label == 1][column].dropna(), s2[s2.label == 0][column].dropna(), **keys,
                                   negatives=neg_type, stratum='complex:at_least_one_side'))
    return rows


def degree_mean(metrics):
    """Mean AUROC/AUPR of the three degree-preserving seeds per (set, model, scorer, region, stratum)."""
    deg = metrics[metrics.negatives.str.startswith('degree')]
    keys = ['set', 'model', 'scorer', 'region', 'stratum']
    out = deg.groupby(keys, dropna=False).agg(auroc=('auroc', 'mean'), aupr=('aupr', 'mean'),
                                              auroc_min_seed=('auroc', 'min'), n_seeds=('auroc', 'size')).reset_index()
    out['negatives'] = 'degree_mean'
    return out


# ---------------------------------------------------------------- ranking
def random_hit(n, m, k):
    """P(best of m true receptors is in the top k of n candidates) under random scores."""
    return 1.0 if k >= n else 1 - comb(n - m, k)/comb(n, k)


def ranking(folder, scores):
    folder = Path(folder)
    pos = pd.read_csv(folder/'positives.csv', dtype=str)
    annotated = pd.read_csv(folder/'pool_annotated.csv', dtype=str)
    acc = dict(zip(pos.ligand, pos.ligand_accessions)) | dict(zip(pos.receptor, pos.receptor_accessions))
    receptors = sorted(set(pos.receptor))
    positive = {frozenset(k) for k in zip(pos.ligand, pos.receptor)}
    other = {frozenset(k) for k, p in zip(zip(annotated.ligand, annotated.receptor), annotated.positive)
             if p != 'True'}
    rows = []
    for l in sorted(set(pos.ligand)):
        cands = [r for r in receptors if r != l and frozenset((l, r)) not in other]
        true = [r for r in cands if frozenset((l, r)) in positive]
        check(true, f'{l}: no true receptor among candidates')
        s = {r: scores.pair(acc[l], acc[r]) for r in cands}
        best = max(s[r] for r in true)
        rank = 1 + sum(1 for r in cands if r not in true and s[r] >= best)
        n, m = len(cands), len(true)
        rows.append(dict(ligand=l, n_candidates=n, n_true=m, best_rank=rank,
                         best_true=max(true, key=lambda r: s[r]), best_true_score=best,
                         random_top1=random_hit(n, m, 1), random_top5=random_hit(n, m, 5)))
    return pd.DataFrame(rows)


def ranking_summary(per_ligand, keys):
    rng = np.random.default_rng(BOOT_SEED)
    r = per_ligand.best_rank.to_numpy()
    stats = dict(top1=lambda x: np.mean(x == 1), top5=lambda x: np.mean(x <= 5), median_rank=np.median)
    row = dict(**keys, ligands=len(r), random_top1=per_ligand.random_top1.mean(),
               random_top5=per_ligand.random_top5.mean(), median_candidates=float(per_ligand.n_candidates.median()))
    boots = [r[rng.integers(len(r), size=len(r))] for _ in range(N_BOOT)]
    for name, f in stats.items():
        row[name] = float(f(r))
        lo, hi = np.percentile([f(b) for b in boots], [2.5, 97.5])
        row[f'{name}_lo'], row[f'{name}_hi'] = float(lo), float(hi)
    return row


# ---------------------------------------------------------------- hub diagnostic
def hub(folder, scores, cellchat_counts):
    folder = Path(folder)
    pos = pd.read_csv(folder/'positives.csv', dtype=str)
    annotated = pd.read_csv(folder/'pool_annotated.csv', dtype=str)
    excluded = {frozenset(k) for k in zip(annotated.ligand, annotated.receptor)}
    acc = dict(zip(pos.ligand, pos.ligand_accessions)) | dict(zip(pos.receptor, pos.receptor_accessions))
    values, seen = {}, set()
    for l in sorted(set(pos.ligand)):
        for r in sorted(set(pos.receptor)):
            key = frozenset((l, r))
            if l == r or key in excluded or key in seen:
                continue
            seen.add(key)
            v = scores.pair(acc[l], acc[r])
            values.setdefault(l, []).append(v)
            values.setdefault(r, []).append(v)
    degree = pd.concat([pos.ligand, pos.receptor]).value_counts()
    rows = [dict(gene=g, mean_unannotated_score=float(np.mean(v)), n_unannotated=len(v),
                 cellchat_rows=int(cellchat_counts.get(g.upper(), 0)), degree_in_set=int(degree.get(g, 0)))
            for g, v in sorted(values.items())]
    return pd.DataFrame(rows)


def spearman_ci(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    rho = float(spearmanr(x, y).statistic)
    rng = np.random.default_rng(BOOT_SEED)
    boots = []
    for _ in range(N_BOOT):
        i = rng.integers(len(x), size=len(x))
        boots.append(spearmanr(x[i], y[i]).statistic)
    lo, hi = np.nanpercentile(boots, [2.5, 97.5])
    return rho, float(lo), float(hi)


def cellchat_gene_counts(data):
    """Rows of the human CellChatDB (all four annotations) in which a gene takes part, as a subunit too."""
    it = pd.read_csv(data/'cellchatdb_human'/'interaction.csv', dtype=str, keep_default_na=False)
    cx = pd.read_csv(data/'cellchatdb_human'/'complex.csv', dtype=str, keep_default_na=False).set_index('complex')
    complexes = {k: [v for v in row if v.strip()] for k, row in cx.iterrows()}
    counts = {}
    for r in it.to_dict('records'):
        genes = set()
        for side, sym in [('ligand', 'ligand.symbol'), ('receptor', 'receptor.symbol')]:
            genes.update(complexes.get(r[side]) or [r[side]])
            genes.update(g.strip() for g in str(r[sym]).split(',') if g.strip())
        for g in {g.upper() for g in genes}:
            counts[g] = counts.get(g, 0) + 1
    return counts


# ---------------------------------------------------------------- consistency
def reverse_id(pid):
    """'A__B' -> 'B__A' (benchmark jobs); 'L__R__LR' -> 'L__R__RL' (diagnostics jobs)."""
    parts = pid.split('__')
    if len(parts) == 3 and parts[2] in ('LR', 'RL'):
        return f"{parts[0]}__{parts[1]}__{'RL' if parts[2] == 'LR' else 'LR'}"
    return f'{parts[1]}__{parts[0]}'


def consistency(official, ours, jobset, model):
    a = official.table.merge(ours.table[['pair_id', 'score']], on='pair_id', suffixes=('_official', '_ours'))
    check(len(a) == len(official.table) == len(ours.table), 'official and own scores cover different jobs')
    a['truncated'] = a.truncated.astype(str) == 'True'
    full = a[~a.truncated]
    diff = (full.score_official - full.score_ours).abs()
    rows = {pid: (o, m) for pid, o, m in zip(full.pair_id, full.score_official, full.score_ours)}
    means = []
    for pid, (o, m) in rows.items():
        rev = reverse_id(pid)
        if pid < rev and rev in rows:
            o2, m2 = rows[rev]
            means.append(abs((o + o2)/2 - (m + m2)/2))
    worst = full.loc[diff.idxmax(), 'pair_id'] if len(full) else None
    return dict(jobset=jobset, model=model, jobs=len(a), truncated_by_official=int(a.truncated.sum()),
                compared=len(full), max_abs_diff_order=float(diff.max()), mean_abs_diff_order=float(diff.mean()),
                max_abs_diff_pair_mean=float(max(means)) if means else np.nan, worst_pair_id=worst,
                within_tolerance=bool(diff.max() <= TOL))


# ---------------------------------------------------------------- training overlap
def seq_key(sequence, cache):
    if sequence not in cache:
        cache[sequence] = hashlib.sha1(sequence.strip().upper().encode()).hexdigest()
    return cache[sequence]


def training_pairs(root, model):
    """Unordered pairs (exact sequences) seen with label 1 / label 0, and every protein seen."""
    seen_pos, seen_neg, proteins, cache = set(), set(), set(), {}
    for name in TRAINING[model]:
        path = Path(root)/name
        for chunk in pd.read_csv(path, usecols=['query', 'text', 'label'], dtype=str, chunksize=200000):
            q = [seq_key(x, cache) for x in chunk['query']]
            t = [seq_key(x, cache) for x in chunk['text']]
            proteins.update(q)
            proteins.update(t)
            for a, b, lab in zip(q, t, chunk.label):
                (seen_pos if lab.strip() == '1' else seen_neg).add(frozenset((a, b)))
    return seen_pos, seen_neg, proteins


def overlap_category(t, sequences, seen_pos, seen_neg, proteins):
    cache = {}
    key = {a: seq_key(s, cache) for a, s in sequences.items()}
    cats = []
    for la, ra in zip(t.ligand_accessions, t.receptor_accessions):
        pairs = [frozenset((key[a], key[b])) for a in la.split(';') for b in ra.split(';') if a != b]
        prots = [key[a] for a in la.split(';') + ra.split(';')]
        if any(p in seen_pos for p in pairs):
            cats.append('seen_positive')
        elif any(p in seen_neg for p in pairs):
            cats.append('seen_negative')
        elif all(p in proteins for p in prots):
            cats.append('proteins_seen')
        else:
            cats.append('protein_unseen')
    return cats


# ---------------------------------------------------------------- figures
def histograms(t, title, path, column='score'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    groups = [('positives', t[t.label == 1][column])]
    groups += [(n, t[t.negative == n][column]) for n in ['random_s1', 'degree_s1', 'degree_s2', 'degree_s3']
               if (t.negative == n).any()]
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8))
    bins = np.linspace(0, 1, 41)
    for name, v in groups:
        axes[0].hist(v.dropna(), bins=bins, histtype='step', density=True, label=f'{name} (n={v.notna().sum()})')
        lv = np.log10(np.clip(v.dropna().to_numpy(float), 1e-8, 1))
        axes[1].hist(lv, bins=np.linspace(-8, 0, 41), histtype='step', density=True, label=name)
    axes[0].set_xlabel('score (mean of both orders)')
    axes[1].set_xlabel('log10 score (clipped at 1e-8)')
    axes[0].set_ylabel('density')
    axes[0].legend(fontsize=7)
    fig.suptitle(title, fontsize=10, x=0.01, ha='left')
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# ---------------------------------------------------------------- evaluate
def evaluate(args):
    runs, data, dataset, out = Path(args.runs), Path(args.data), Path(args.dataset), Path(args.output)
    check(not out.exists(), f'{out} exists')
    out.mkdir(parents=True)
    (out/'pair_scores').mkdir()
    (out/'figures').mkdir()
    load = lambda js, m, s: Scores(runs/'scores'/f'{js}__{m}__{s}.csv.gz')
    metrics, rank_rows, hub_rows, cons_rows = [], [], [], []
    cc_counts = cellchat_gene_counts(data)
    sequences = {}
    for species in ['human', 'mouse']:
        for p in (data/f'proteins_{species}'/'uniprot').glob('*.json'):
            sequences[p.stem] = json.loads(p.read_text())['sequence']['value']
    overlap = {}
    for model in MODELS:
        try:
            overlap[model] = training_pairs(data/'plm_training', model)
        except (FileNotFoundError, ValueError) as exc:
            overlap[model] = str(exc)
    for set_name, (jobset, models) in SCOPE.items():
        for model in models:
            ours, official = load(jobset, model, 'ours'), load(jobset, model, 'official')
            extra = dict(score_official=official)
            if set_name == 'H-M' and model == MAIN_MODEL:
                extra['score_ecd'] = load('human_ecd', model, 'ours')
            t = score_table(dataset/set_name, ours, extra)
            keys = dict(set=set_name, model=model, scorer='ours', region='full')
            metrics += set_metrics(t, keys)
            if set_name == 'H-A':
                metrics += [r for r in set_metrics(t, dict(keys, scorer='official'), 'score_official')
                            if r['stratum'].startswith('length:other')]
            if 'score_ecd' in t:
                both = t[t.score_ecd.notna()]
                metrics += [r for r in set_metrics(both, dict(keys, region='full_on_ecd_subset')) if r['stratum'] == 'all']
                metrics += [r for r in set_metrics(both, dict(keys, region='extracellular'), 'score_ecd')
                            if r['stratum'] == 'all']
            if isinstance(overlap[model], tuple):
                t['training'] = overlap_category(t, sequences, *overlap[model])
                for neg_type in sorted(n for n in t.negative.unique() if n != 'positive'):
                    neg = t[t.negative == neg_type]
                    for cat in ['seen_positive', 'seen_negative', 'proteins_seen', 'protein_unseen']:
                        p = t[(t.label == 1) & (t.training == cat)].score
                        metrics.append(metric_row(p, neg.score, **keys, negatives=neg_type, stratum=f'training:{cat}'))
                    unseen_pos = t[(t.label == 1) & (t.training != 'seen_positive')].score
                    unseen_neg = neg[~neg.training.isin(['seen_positive', 'seen_negative'])].score
                    metrics.append(metric_row(unseen_pos, unseen_neg, **keys, negatives=neg_type,
                                              stratum='training:pos_not_seen_vs_neg_not_seen'))
            t.drop(columns=['query', 'text'], errors='ignore').to_csv(
                out/'pair_scores'/f'{set_name}__{model}.csv.gz', index=False)
            histograms(t, f'{set_name}, {model}, plm_score.py, full length', out/'figures'/f'hist_{set_name}_{model}.png')
            if set_name.endswith('-M'):
                per = ranking(dataset/set_name, ours)
                per.to_csv(out/'pair_scores'/f'ranking_{set_name}__{model}.csv', index=False)
                rank_rows.append(ranking_summary(per, dict(set=set_name, model=model, scorer='ours')))
            if set_name == 'H-M':
                h = hub(dataset/set_name, ours, cc_counts)
                h.to_csv(out/'pair_scores'/f'hub_{set_name}__{model}.csv', index=False)
                for other in ['cellchat_rows', 'degree_in_set']:
                    rho, lo, hi = spearman_ci(h.mean_unannotated_score, h[other])
                    hub_rows.append(dict(set=set_name, model=model, versus=other, proteins=len(h),
                                         spearman=rho, lo=lo, hi=hi))
        for model in models:
            if not any(r['jobset'] == jobset and r['model'] == model for r in cons_rows):
                cons_rows.append(consistency(load(jobset, model, 'official'), load(jobset, model, 'ours'),
                                             jobset, model))
    m = pd.DataFrame(metrics)
    m = pd.concat([m, degree_mean(m)], ignore_index=True)
    m.to_csv(out/'metrics.csv', index=False)
    ranks = pd.DataFrame(rank_rows)
    ranks.to_csv(out/'ranking.csv', index=False)
    pd.DataFrame(hub_rows).to_csv(out/'hub.csv', index=False)
    pd.DataFrame(cons_rows).to_csv(out/'consistency.csv', index=False)
    main = m[(m.set == MAIN_SET) & (m.model == MAIN_MODEL) & (m.scorer == 'ours') & (m.region == 'full')
             & (m.stratum == 'all')].set_index('negatives')
    top5 = float(ranks[(ranks.set == MAIN_SET) & (ranks.model == MAIN_MODEL)].top5.iloc[0])
    checks = dict(random_auroc=float(main.loc['random_s1', 'auroc']),
                  degree_auroc=float(main.loc['degree_mean', 'auroc']), top5=top5)
    verdict = dict(criteria=CRITERIA, values=checks,
                   met={k: checks[k] >= CRITERIA[k] for k in CRITERIA},
                   degree_each_seed={s: float(main.loc[s, 'auroc']) for s in ['degree_s1', 'degree_s2', 'degree_s3']},
                   training_overlap_unavailable={k: v for k, v in overlap.items() if not isinstance(v, tuple)})
    verdict['usable'] = all(verdict['met'].values())
    write_json(out/'verdict.json', dict(status='complete', created_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'),
                                        main=dict(set=MAIN_SET, model=MAIN_MODEL, scorer='plm_score.py', region='full'),
                                        **verdict))
    print(json.dumps(verdict, indent=1))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    c = sub.add_parser('collect')
    c.add_argument('--runs', default='plm_lr_benchmark_runs')
    e = sub.add_parser('evaluate')
    e.add_argument('--runs', default='plm_lr_benchmark_runs')
    e.add_argument('--data', default='data/lr_benchmark')
    e.add_argument('--dataset', default='reports/plm_lr_benchmark/dataset')
    e.add_argument('--output', default='reports/plm_lr_benchmark/results')
    return p


if __name__ == '__main__':
    a = parser().parse_args()
    {'collect': collect, 'evaluate': evaluate}[a.command](a)
