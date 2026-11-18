"""Isoform / paralog ranking benchmark for PLM-interact scores.

Plan: reports/isoform_lr_bench_01/plan.md (commit it before the main run).
Environment: commot-research (numpy, pandas). Scoring is done by plm_score.py in the
plminteract environment; this file never imports torch.

  jobs      dataset -> <out>/plm_input.csv (pair_id,query,text; both orientations) + index.csv
  merge     plm_score.py output + index -> one score per dataset row (mean of the two orientations)
  evaluate  scored rows -> metrics.json, per_contrast.csv, pairwise.csv, row_scores.csv, verdict.json
  selftest  synthetic scores (oracle / random / canonical-only / longest-only) through evaluate;
            checks that the pre-registered decision rule behaves as intended. CPU only.

Every path is a command-line argument; nothing is written unless an output path is given, and an
existing output is never overwritten.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

# ------------------------------------------------------------------ pre-registered constants
DATASET_SHA256 = 'd3738f564f48fa5792e7e0ee0d2e7f2e5246407a83e754a975cdca3f84bedb3b'
EXPECTED_CHECKPOINTS = {
    'humanV12': dict(checkpoint_name='danliu1226/PLM-interact-650M-humanV12',
                     checkpoint_sha256='92565937d65a9b91a4e897c99653c972b99cfe37067eb42d9e783af688350547'),
}
PRIMARY_MODEL = 'humanV12'
RANKED_ROLES = ('reference', 'competitor')
TIERS = ('H_isoform', 'E_paralog')
IDENTITY_BINS = [0, 50, 70, 90, 99, 100.001]
OFFICIAL_MAX_TOKENS = 1603
CRITERIA = dict(alpha=0.05, H_min_hit_at_1=0.45, H_max_mean_norm_rank=0.40,
                H_min_margin_over_trivial=0.10, E_min_hit_at_1=0.50, unseen_min_contrasts=30)


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False, default=str))


def fresh_dir(path):
    path = Path(path)
    if path.exists() and any(path.iterdir()):
        raise FileExistsError(f'{path} exists and is not empty; use a new output directory')
    path.mkdir(parents=True, exist_ok=True)
    return path


def load_dataset(path, check_sha=True):
    sha = digest(path)
    if check_sha and sha != DATASET_SHA256:
        raise ValueError(f'{path}: sha256 {sha} != pre-registered {DATASET_SHA256}. '
                         'Pass --exploratory to run on a different file (results are then exploratory).')
    data = pd.read_csv(path, dtype={'pair_id': str})
    if not data.pair_id.is_unique:
        raise ValueError('pair_id must be unique')
    return data, sha


# ------------------------------------------------------------------ jobs
def select_contrasts(data, limit, seed, include_longest):
    if not limit:
        return data
    rng = np.random.default_rng(seed)
    keep = []
    for tier in TIERS:
        ids = sorted(data.loc[data.tier == tier, 'contrast_id'].unique())
        keep += list(rng.choice(ids, size=min(limit, len(ids)), replace=False))
    if include_longest:
        keep.append(data.loc[data.n_tokens.idxmax(), 'contrast_id'])
    return data[data.contrast_id.isin(set(keep))]


def jobs(args):
    data, sha = load_dataset(args.dataset, not args.exploratory)
    data = select_contrasts(data, args.limit_contrasts, args.seed, args.include_longest)
    out = fresh_dir(args.output_dir)
    rows, index = [], []
    for r in data.itertuples(index=False):
        for orient, q, t in (('vp', r.seq_variant, r.seq_partner), ('pv', r.seq_partner, r.seq_variant)):
            job = f'{r.pair_id}__{orient}'
            rows.append(dict(pair_id=job, query=q, text=t))
            index.append(dict(job_id=job, pair_id=r.pair_id, orientation=orient))
    pd.DataFrame(rows).to_csv(out/'plm_input.csv', index=False)
    pd.DataFrame(index).to_csv(out/'index.csv', index=False)
    info = dict(created_utc=now(), dataset=str(Path(args.dataset).resolve()), dataset_sha256=sha,
                exploratory=bool(args.exploratory), limit_contrasts=args.limit_contrasts, seed=args.seed,
                include_longest=bool(args.include_longest), n_rows=len(data),
                n_contrasts=int(data.contrast_id.nunique()), n_jobs=len(rows),
                n_rows_over_official_max_tokens=int((data.n_tokens > OFFICIAL_MAX_TOKENS).sum()),
                max_tokens=int(data.n_tokens.max()),
                plm_input_sha256=digest(out/'plm_input.csv'))
    write_json(out/'jobs.json', info)
    print(json.dumps(info, indent=2))


# ------------------------------------------------------------------ merge
def merge(args):
    data, sha = load_dataset(args.dataset, not args.exploratory)
    jobs_dir = Path(args.jobs_dir)
    index = pd.read_csv(jobs_dir/'index.csv', dtype=str)
    jobs_info = json.loads((jobs_dir/'jobs.json').read_text())
    if jobs_info['dataset_sha256'] != sha:
        raise ValueError('jobs were built from a different dataset file')
    meta = json.loads(Path(str(args.scores) + '.meta.json').read_text())
    if meta.get('status') != 'complete':
        raise ValueError('score file is not marked complete')
    if meta['input_sha256'] != digest(jobs_dir/'plm_input.csv'):
        raise ValueError('scores were produced from a different plm_input.csv')
    scores = pd.read_csv(args.scores, dtype={'pair_id': str})
    if not scores.pair_id.is_unique:
        raise ValueError('a job was scored more than once')
    missing = set(index.job_id) - set(scores.pair_id)
    extra = set(scores.pair_id) - set(index.job_id)
    if missing or extra:
        raise ValueError(f'{len(missing)} jobs unscored, {len(extra)} unexpected scores')
    s = index.merge(scores.rename(columns={'pair_id': 'job_id'}), on='job_id')
    wide = s.pivot(index='pair_id', columns='orientation', values='score')
    wide.columns = [f'score_{c}' for c in wide.columns]
    wide['score'] = wide[['score_vp', 'score_pv']].mean(axis=1)
    out = data.merge(wide.reset_index(), on='pair_id', how='inner')
    if len(out) != index.pair_id.nunique():
        raise ValueError('row count mismatch after merge')
    path = Path(args.output)
    if path.exists():
        raise FileExistsError(f'{path} exists')
    path.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(path, index=False)
    keep = ['checkpoint_name', 'checkpoint_revision', 'checkpoint_sha256', 'base_model', 'torch',
            'transformers', 'device', 'started_utc', 'finished_utc', 'max_tokens', 'scored_beyond_max_tokens']
    write_json(str(path) + '.meta.json', dict(
        created_utc=now(), dataset_sha256=sha, exploratory=bool(args.exploratory or jobs_info['exploratory']),
        jobs=jobs_info, scores=str(Path(args.scores).resolve()), scores_sha256=digest(args.scores),
        n_rows=len(out), output_sha256=digest(path),
        **{k: (len(meta[k]) if k == 'scored_beyond_max_tokens' else meta.get(k)) for k in keep}))
    print(f'Merged {len(out)} rows ({len(index)} jobs) -> {path}')


# ------------------------------------------------------------------ metrics
def ranks_desc(values):
    """1 = highest. Ties share the average rank."""
    v = np.asarray(values, float)
    order = np.argsort(-v, kind='mergesort')
    r = np.empty(len(v))
    r[order] = np.arange(1, len(v) + 1)
    for val in np.unique(v):
        m = v == val
        if m.sum() > 1:
            r[m] = r[m].mean()
    return r


def contrast_table(rows, score_col='score'):
    out = []
    for cid, v in rows.groupby('contrast_id', sort=True):
        ref = (v.role_in_contrast == 'reference').to_numpy()
        if ref.sum() != 1 or ref.all():
            raise ValueError(f'{cid}: needs one reference and >=1 competitor')
        s = v[score_col].to_numpy(float)
        k = len(v)
        rank = ranks_desc(s)[ref][0]
        comp = v[~ref]
        out.append(dict(contrast_id=cid, tier=v.tier.iloc[0], reference_id=v.variant_id[ref].iloc[0],
                        reference_gene=v.variant_gene[ref].iloc[0], partner_id=v.partner_id.iloc[0], k=k,
                        rank=rank, norm_rank=(rank - 1) / (k - 1), hit_at_1=int(rank == 1),
                        score_reference=float(s[ref][0]), score_best_competitor=float(s[~ref].max()),
                        delta=float(s[ref][0] - s[~ref].max()),
                        min_identity=float(comp.pct_id_to_reference.min()),
                        max_tokens=int(v.n_tokens.max())))
    return pd.DataFrame(out)


def permutation_p(table, n_perm, rng):
    """Null: the reference takes each of the k ranks with equal probability, independently per contrast."""
    k = table.k.to_numpy()
    sims = rng.integers(0, k, size=(n_perm, len(k)))          # 0 = top
    null_nr = (sims / (k - 1)).mean(axis=1)
    null_h1 = (sims == 0).mean(axis=1)
    obs_nr, obs_h1 = table.norm_rank.mean(), table.hit_at_1.mean()
    return (float((1 + (null_nr <= obs_nr).sum()) / (1 + n_perm)),
            float((1 + (null_h1 >= obs_h1).sum()) / (1 + n_perm)))


def cluster_ci(table, col, n_boot, rng, cluster='reference_gene'):
    groups = [g[col].to_numpy(float) for _, g in table.groupby(cluster)]
    if len(groups) < 2 or n_boot <= 0:
        return [None, None]
    stats = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        stats.append(np.concatenate([groups[i] for i in pick]).mean())
    return [round(float(np.quantile(stats, 0.025)), 4), round(float(np.quantile(stats, 0.975)), 4)]


def summarise(table, n_perm, n_boot, rng):
    if not len(table):
        return dict(n_contrasts=0)
    p_nr, p_h1 = permutation_p(table, n_perm, rng)
    return dict(n_contrasts=int(len(table)), n_reference_genes=int(table.reference_gene.nunique()),
                hit_at_1=round(float(table.hit_at_1.mean()), 4),
                hit_at_1_ci95_gene_bootstrap=cluster_ci(table, 'hit_at_1', n_boot, rng),
                chance_hit_at_1=round(float((1 / table.k).mean()), 4),
                mean_norm_rank=round(float(table.norm_rank.mean()), 4),
                mean_norm_rank_ci95_gene_bootstrap=cluster_ci(table, 'norm_rank', n_boot, rng),
                p_hit_at_1=round(p_h1, 5), p_norm_rank=round(p_nr, 5),
                median_delta=round(float(table.delta.median()), 4),
                frac_delta_gt_0=round(float((table.delta > 0).mean()), 4))


def pairwise_table(rows, score_col='score'):
    out = []
    for cid, v in rows.groupby('contrast_id', sort=True):
        ref = v[v.role_in_contrast == 'reference'].iloc[0]
        for c in v[v.role_in_contrast == 'competitor'].itertuples():
            sc = getattr(c, score_col)
            out.append(dict(contrast_id=cid, tier=c.tier, reference_id=ref.variant_id, competitor_id=c.variant_id,
                            pct_id_to_reference=c.pct_id_to_reference, score_reference=ref[score_col],
                            score_competitor=sc, win=1.0 if ref[score_col] > sc else (0.5 if ref[score_col] == sc else 0.0)))
    return pd.DataFrame(out)


def by_identity(pw):
    out = {}
    cut = pd.cut(pw.pct_id_to_reference, IDENTITY_BINS, right=False)
    for b, g in pw.groupby(cut, observed=False):
        out[f'{b.left:g}-{min(b.right, 100):g}'] = dict(n_pairs=int(len(g)),
            win_rate=round(float(g.win.mean()), 4) if len(g) else None,
            n_contrasts=int(g.contrast_id.nunique()))
    return out


def trivial_baselines(rows):
    h = rows[rows.tier == 'H_isoform'].copy()
    h['canonical'] = (h.variant_status == 'Displayed').astype(float)
    h['longest'] = h.len_variant.astype(float)
    return {name: round(float(contrast_table(h, name).hit_at_1.mean()), 4) for name in ('canonical', 'longest')}


def training_overlap(rows, training_root, model):
    """Same rule as lr_eval.py: exact-sequence SHA1, unordered pair seen with label 1."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from lr_eval import training_pairs, seq_key
    seen_pos, seen_neg, proteins = training_pairs(training_root, model)
    cache = {}
    pair = [frozenset((seq_key(a, cache), seq_key(b, cache))) for a, b in zip(rows.seq_variant, rows.seq_partner)]
    rows = rows.copy()
    rows['train_pos'] = [p in seen_pos for p in pair]
    rows['train_neg'] = [p in seen_neg for p in pair]
    rows['variant_seen'] = [seq_key(s, cache) in proteins for s in rows.seq_variant]
    return rows


def decide(H, E, H_unseen, baselines):
    c = CRITERIA
    crit = {}
    crit['C1_isoform_above_chance'] = bool(H['p_hit_at_1'] < c['alpha'] and H['p_norm_rank'] < c['alpha'])
    crit['C2_isoform_effect_size'] = bool(H['hit_at_1'] >= c['H_min_hit_at_1'] and
                                          H['mean_norm_rank'] <= c['H_max_mean_norm_rank'])
    best_trivial = max(baselines.values())
    crit['C3_beats_trivial_baselines'] = bool(H['hit_at_1'] - best_trivial >= c['H_min_margin_over_trivial'])
    if H_unseen is None:
        crit['C5_unseen_subset'] = 'not_run'
    elif H_unseen.get('n_contrasts', 0) < c['unseen_min_contrasts']:
        crit['C5_unseen_subset'] = 'not_evaluable'
    else:
        crit['C5_unseen_subset'] = bool(H_unseen['p_norm_rank'] < c['alpha'])
    crit['C4_paralog_tier'] = bool(E['hit_at_1'] >= c['E_min_hit_at_1'] and E['p_hit_at_1'] < c['alpha'])
    c5 = crit['C5_unseen_subset']
    core = crit['C1_isoform_above_chance'] and crit['C2_isoform_effect_size'] and crit['C3_beats_trivial_baselines']
    if c5 == 'not_run':
        iso = 'incomplete (training overlap not run)'
    elif c5 is False:
        iso = 'FAIL'
    else:
        iso = 'PASS' if core else 'FAIL'
    crit['VERDICT_isoform_resolution'] = iso
    crit['VERDICT_paralog_ranking'] = 'PASS' if crit['C4_paralog_tier'] else 'FAIL'
    crit['best_trivial_hit_at_1'] = best_trivial
    return crit


def evaluate_frame(scored, n_perm, n_boot, seed, training_root=None, model=PRIMARY_MODEL):
    rng = np.random.default_rng(seed)
    rows = scored[scored.role_in_contrast.isin(RANKED_ROLES)].copy()
    if training_root:
        rows = training_overlap(rows, training_root, model)
    tables, summary = {}, {}
    for tier in TIERS:
        tables[tier] = contrast_table(rows[rows.tier == tier])
        summary[tier] = summarise(tables[tier], n_perm, n_boot, rng)
    per_contrast = pd.concat(tables.values(), ignore_index=True)
    H_unseen = None
    if training_root:
        ref_seen = rows[rows.role_in_contrast == 'reference'].set_index('contrast_id').train_pos
        per_contrast['reference_pair_in_training_pos'] = per_contrast.contrast_id.map(ref_seen)
        any_comp = rows[rows.role_in_contrast == 'competitor'].groupby('contrast_id').train_pos.any()
        per_contrast['any_competitor_pair_in_training_pos'] = per_contrast.contrast_id.map(any_comp)
        h = per_contrast[per_contrast.tier == 'H_isoform']
        H_unseen = summarise(h[~h.reference_pair_in_training_pos.astype(bool)], n_perm, n_boot, rng)
        summary['H_isoform_reference_not_in_training'] = H_unseen
        summary['training_overlap'] = {t: dict(
            contrasts_reference_pair_in_training_pos=int(per_contrast[per_contrast.tier == t].reference_pair_in_training_pos.sum()),
            contrasts_any_competitor_pair_in_training_pos=int(per_contrast[per_contrast.tier == t].any_competitor_pair_in_training_pos.sum()),
            ranked_rows_variant_seen=int(rows[rows.tier == t].variant_seen.sum()),
            ranked_rows=int((rows.tier == t).sum())) for t in TIERS}
    pw = pairwise_table(rows)
    summary['H_pairwise_by_identity'] = by_identity(pw[pw.tier == 'H_isoform'])
    summary['pairwise_win_rate'] = {t: round(float(pw[pw.tier == t].win.mean()), 4) for t in TIERS}
    h = per_contrast[per_contrast.tier == 'H_isoform']
    summary['H_by_tokens_descriptive'] = {
        f'max_tokens<={OFFICIAL_MAX_TOKENS}': summarise(h[h.max_tokens <= OFFICIAL_MAX_TOKENS], 2000, 0, rng),
        f'max_tokens>{OFFICIAL_MAX_TOKENS}': summarise(h[h.max_tokens > OFFICIAL_MAX_TOKENS], 2000, 0, rng)}
    baselines = trivial_baselines(rows)
    summary['trivial_baselines_H_hit_at_1'] = baselines
    summary['preregistered_criteria'] = decide(summary['H_isoform'], summary['E_paralog'], H_unseen, baselines)
    return summary, per_contrast, pw, rows


def evaluate(args):
    scored = pd.read_csv(args.scored, dtype={'pair_id': str})
    meta = json.loads(Path(str(args.scored) + '.meta.json').read_text())
    exploratory = bool(args.exploratory or meta.get('exploratory'))
    problems = []
    if meta.get('dataset_sha256') != DATASET_SHA256:
        problems.append('dataset differs from the pre-registered file')
    exp = EXPECTED_CHECKPOINTS.get(args.model)
    if exp is None or meta.get('checkpoint_sha256') != exp['checkpoint_sha256']:
        problems.append(f'checkpoint is not the pre-registered {args.model} file')
    if not args.training_root:
        problems.append('training overlap (C5) not run')
    if problems and not exploratory:
        raise ValueError('; '.join(problems) + '. Pass --exploratory to proceed (results are then exploratory).')
    out = fresh_dir(args.output_dir)
    summary, per_contrast, pw, rows = evaluate_frame(scored, args.n_perm, args.n_boot, args.seed,
                                                     args.training_root, args.model)
    head = dict(created_utc=now(), status='exploratory' if exploratory else 'pre-registered',
                problems=problems, scored=str(Path(args.scored).resolve()), scored_sha256=digest(args.scored),
                model=args.model, checkpoint_name=meta.get('checkpoint_name'),
                checkpoint_sha256=meta.get('checkpoint_sha256'), dataset_sha256=meta.get('dataset_sha256'),
                n_perm=args.n_perm, n_boot=args.n_boot, seed=args.seed, criteria=CRITERIA,
                training_root=str(Path(args.training_root).resolve()) if args.training_root else None)
    write_json(out/'metrics.json', dict(**head, **summary))
    write_json(out/'verdict.json', dict(**head, **summary['preregistered_criteria']))
    per_contrast.to_csv(out/'per_contrast.csv', index=False)
    pw.to_csv(out/'pairwise.csv', index=False)
    drop = [c for c in ('seq_variant', 'seq_partner') if c in rows]
    scored.drop(columns=drop).to_csv(out/'row_scores.csv', index=False)
    print(json.dumps(dict(status=head['status'], problems=problems, H=summary['H_isoform'], E=summary['E_paralog'],
                          trivial=summary['trivial_baselines_H_hit_at_1'],
                          criteria=summary['preregistered_criteria']), indent=2))


# ------------------------------------------------------------------ selftest
def synthetic(data, kind, rng):
    d = data.copy()
    ref = (d.role_in_contrast == 'reference').to_numpy()
    if kind == 'oracle':
        d['score'] = np.where(ref, 0.9, rng.uniform(0, 0.5, len(d)))
    elif kind == 'random':
        d['score'] = rng.uniform(0, 1, len(d))
    elif kind == 'canonical_only':
        d['score'] = np.where(d.variant_status == 'Displayed', 0.8, 0.2) + rng.uniform(0, 1e-3, len(d))
        e = (d.tier == 'E_paralog').to_numpy()
        d.loc[e, 'score'] = np.where(ref[e], 0.9, rng.uniform(0, 0.5, e.sum()))   # gene-level signal present
    elif kind == 'longest_only':
        d['score'] = d.len_variant / d.len_variant.max()
    return d


SELFTEST_EXPECT = {
    'oracle': dict(VERDICT_isoform_resolution='PASS', VERDICT_paralog_ranking='PASS'),
    'random': dict(VERDICT_isoform_resolution='FAIL', VERDICT_paralog_ranking='FAIL'),
    'canonical_only': dict(VERDICT_isoform_resolution='FAIL', VERDICT_paralog_ranking='PASS',
                           C3_beats_trivial_baselines=False),
    'longest_only': dict(C3_beats_trivial_baselines=False),
}


def selftest(args):
    data, sha = load_dataset(args.dataset, not args.exploratory)
    out = fresh_dir(args.output_dir)
    rng = np.random.default_rng(args.seed)
    report, ok = {}, True
    for kind, expect in SELFTEST_EXPECT.items():
        summary, *_ = evaluate_frame(synthetic(data, kind, rng), args.n_perm, 200, args.seed)
        crit = summary['preregistered_criteria']
        # no training root here: C5 is 'not_run', so judge the core verdict
        core = 'PASS' if (crit['C1_isoform_above_chance'] and crit['C2_isoform_effect_size']
                          and crit['C3_beats_trivial_baselines']) else 'FAIL'
        got = dict(crit, VERDICT_isoform_resolution=core)
        checks = {k: got[k] == v for k, v in expect.items()}
        ok &= all(checks.values())
        report[kind] = dict(H_hit_at_1=summary['H_isoform']['hit_at_1'], E_hit_at_1=summary['E_paralog']['hit_at_1'],
                            criteria=got, expected=expect, passed=all(checks.values()))
    write_json(out/'selftest.json', dict(created_utc=now(), dataset_sha256=sha, passed=bool(ok), cases=report))
    print(json.dumps({k: dict(passed=v['passed'], H=v['H_hit_at_1'], E=v['E_hit_at_1']) for k, v in report.items()},
                     indent=2))
    if not ok:
        raise SystemExit('selftest FAILED')
    print('selftest passed')


# ------------------------------------------------------------------ cli
def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    j = sub.add_parser('jobs')
    j.add_argument('--dataset', required=True); j.add_argument('--output-dir', required=True)
    j.add_argument('--limit-contrasts', type=int, default=0, help='per tier; 0 = all (smoke test: 3)')
    j.add_argument('--include-longest', action='store_true', help='also take the contrast with the longest pair')
    j.add_argument('--seed', type=int, default=0)
    m = sub.add_parser('merge')
    m.add_argument('--dataset', required=True); m.add_argument('--jobs-dir', required=True)
    m.add_argument('--scores', required=True, help='plm_score.py output (with .meta.json)')
    m.add_argument('--output', required=True)
    e = sub.add_parser('evaluate')
    e.add_argument('--scored', required=True, help='merge output (with .meta.json)')
    e.add_argument('--output-dir', required=True)
    e.add_argument('--training-root', help='data/lr_benchmark/plm_training (as in lr_eval.py)')
    e.add_argument('--model', default=PRIMARY_MODEL, choices=['humanV11', 'humanV12', 'LeakageFree'])
    e.add_argument('--n-perm', type=int, default=10000); e.add_argument('--n-boot', type=int, default=2000)
    e.add_argument('--seed', type=int, default=20261118)
    s = sub.add_parser('selftest')
    s.add_argument('--dataset', required=True); s.add_argument('--output-dir', required=True)
    s.add_argument('--n-perm', type=int, default=2000); s.add_argument('--seed', type=int, default=7)
    for x in (j, m, e, s):
        x.add_argument('--exploratory', action='store_true',
                       help='allow a non-registered dataset/checkpoint or a missing check; output marked exploratory')
    return p


if __name__ == '__main__':
    a = parser().parse_args()
    dict(jobs=jobs, merge=merge, evaluate=evaluate, selftest=selftest)[a.command](a)
