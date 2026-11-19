"""Experiment 1 (reports/plm_length_01/plan.md): does humanV12 separate long pairs worse than short pairs?

  fetch     download the Bernett et al. (2024) test split at a pinned HF revision and check its SHA256
  prepare   test split + V12 training set -> primary set (pairs not seen in V12 training), length groups,
            training-overlap strata, the E2 sample with its padded sequences, scoring jobs in chunks
  work      worker w of n: score every n-th chunk with plm_score.py; finished chunks are skipped (resumable)
  merge     chunk scores -> one score per unordered pair and variant (mean of the two input orders)
  evaluate  E1 / E2 statistics and the pre-registered decisions (plan sections 4-6); LR replication (descriptive)

Environments: commot-research for fetch / prepare / merge / evaluate; the PLM-interact environment for scoring,
reached through --plm-python (work). All paths are command-line arguments.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent
BERNETT_REPO, BERNETT_REV = 'danliu1226/Bernett_benchmarking', '5d2ad03baa165c27df32a2eadf066462a2a83073'
TEST_FILE = 'pairs_uniprot_seqs_test.csv'
TEST_SHA256 = '89f164f9bd391fd3a2cb400cb3b402d87b73af8142143c58cac56d55d3d7abcc'
TRAIN_REL = 'STRING_V12_TrainingSet/protein.pairs_human_V12_train.tsv'            # under --training-root
TRAIN_SHA256 = 'fc3c0b5df94ee0cd0ae1934fdefb49f2df72c21e83ea2244aca2ad6af8376f43'  # @819f938
CHECKPOINT = ('danliu1226/PLM-interact-650M-humanV12', '2d0e7fb3deef94efa3f217daccc79092a6a42ceb')
ESM_LETTERS = set('ACDEFGHIKLMNPQRSTVWYXBUZO')
MAX_TOKENS = 1603                      # official training max_length
BINS = [0, 600, 1000, 1603, 2400, np.inf]
BIN_LABELS = ['<=600', '601-1000', '1001-1603', '1604-2400', '>2400']
STRATA = ['both_seen', 'one_seen', 'none_seen']
E2_N, E2_MAX_TOKENS, TARGETS, PADS = 2000, 1000, (1603, 2400, 3200), ('shuffle', 'foreign')
SEED, N_BOOT, N_SPLITS, MARGIN, CHUNK = 1, 2000, 20, 0.05, 2000


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def check(cond, msg):
    if not cond:
        raise ValueError(msg)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, default=lambda x: x.item() if hasattr(x, 'item') else str(x)))


def fresh(path):
    path = Path(path)
    check(not path.exists(), f'{path} exists; use a new output')
    path.mkdir(parents=True)
    return path


def seq_key(sequence):
    """Same key as lr_eval.seq_key (exact sequence match)."""
    return hashlib.sha1(sequence.strip().upper().encode()).hexdigest()


# ------------------------------------------------------------------ fetch
def fetch(args):
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    out.parent.mkdir(parents=True, exist_ok=True)
    url = f'https://huggingface.co/datasets/{BERNETT_REPO}/resolve/{BERNETT_REV}/{TEST_FILE}'
    tmp = out.with_suffix('.part')
    urllib.request.urlretrieve(url, tmp)
    sha = digest(tmp)
    check(sha == TEST_SHA256, f'SHA256 {sha} differs from the pinned {TEST_SHA256}')
    tmp.rename(out)
    write_json(Path(str(out) + '.json'), dict(url=url, revision=BERNETT_REV, sha256=sha, fetched_utc=now()))
    print(f'{out}: sha256 ok')


# ------------------------------------------------------------------ prepare
def training_sets(path):
    """Unordered pairs seen with label 1 / 0, every protein seen, and the longest training protein (aa)."""
    seen_pos, seen_neg, proteins, longest = set(), set(), set(), 0
    for chunk in pd.read_csv(path, usecols=['query', 'text', 'label'], dtype=str, chunksize=200000):
        q = [seq_key(x) for x in chunk['query']]
        t = [seq_key(x) for x in chunk['text']]
        longest = max(longest, int(chunk['query'].str.len().max()), int(chunk['text'].str.len().max()))
        proteins.update(q); proteins.update(t)
        for a, b, lab in zip(q, t, chunk.label):
            (seen_pos if lab.strip() == '1' else seen_neg).add(frozenset((a, b)))
    return seen_pos, seen_neg, proteins, longest


def overlap(ka, kb, seen_pos, seen_neg, proteins):
    if frozenset((ka, kb)) in seen_pos or frozenset((ka, kb)) in seen_neg:
        return 'pair_seen'
    return {2: 'both_seen', 1: 'one_seen', 0: 'none_seen'}[(ka in proteins) + (kb in proteins)]


def counts(df, cols):
    """Group sizes with JSON-safe keys, e.g. {'long|1': 3909}."""
    return {'|'.join(map(str, k if isinstance(k, tuple) else (k,))): int(v) for k, v in df.groupby(cols).size().items()}


def allocate(counts, total):
    """Proportional allocation with largest remainders (deterministic)."""
    raw = {k: total * v / sum(counts.values()) for k, v in counts.items()}
    out = {k: int(np.floor(x)) for k, x in raw.items()}
    for k in sorted(raw, key=lambda k: (-(raw[k] - out[k]), k))[:total - sum(out.values())]:
        out[k] += 1
    return out


def e2_sample(primary, n, rng):
    pool = primary[primary.tokens <= E2_MAX_TOKENS]
    picked = []
    for label in (1, 0):
        sub = pool[pool.label == label]
        alloc = allocate(sub.overlap.value_counts().to_dict(), min(n, len(sub)))
        for stratum in sorted(alloc):
            ids = np.sort(sub.uid[sub.overlap == stratum].to_numpy())
            picked += list(rng.choice(ids, size=alloc[stratum], replace=False))
    return set(picked)


def pad_streams(row, seqs, partners, donors, pad, need_a, need_b, number):
    """Residues appended to A and B (lengths need_a, need_b); deterministic per pair and pad type."""
    rng = np.random.default_rng([SEED, PADS.index(pad), number])
    if pad == 'shuffle':
        a, b = seqs[row.Uniprot_a], seqs[row.Uniprot_b]
        return (''.join(rng.choice(list(a), size=need_a, replace=True)),
                ''.join(rng.choice(list(b), size=need_b, replace=True)))
    banned = {row.Uniprot_a, row.Uniprot_b} | partners.get(row.Uniprot_a, set()) | partners.get(row.Uniprot_b, set())
    allowed = [d for d in donors if d not in banned]
    order = rng.permutation(len(allowed))
    stream, k = '', 0
    first = seqs[allowed[order[0]]]
    offset = int(rng.integers(len(first)))
    stream = first[offset:]
    while len(stream) < need_a + need_b:
        k += 1
        stream += seqs[allowed[order[k % len(order)]]]
    return stream[:need_a], stream[need_a:need_a + need_b]


def prepare(args):
    out = fresh(args.output_dir)
    test, train = Path(args.test_file), Path(args.training_root)/TRAIN_REL
    if not args.skip_sha:
        check(digest(test) == TEST_SHA256, f'{test}: SHA256 differs from the pinned file')
        check(digest(train) == TRAIN_SHA256, f'{train}: SHA256 differs from the pinned file')
    data = pd.read_csv(test, dtype={'Uniprot_a': str, 'Uniprot_b': str, 'query': str, 'text': str})
    check({'Uniprot_a', 'Uniprot_b', 'label', 'query', 'text'} <= set(data), 'Unexpected test-file columns')
    data['query'] = data['query'].str.strip().str.upper(); data['text'] = data['text'].str.strip().str.upper()
    ok = data['query'].map(lambda s: set(s) <= ESM_LETTERS) & data['text'].map(lambda s: set(s) <= ESM_LETTERS)
    check(ok.all(), f'{(~ok).sum()} pairs contain letters outside the ESM vocabulary')
    data['uid'] = [f'B{i:05d}' for i in range(len(data))]
    seen_pos, seen_neg, proteins, longest = training_sets(train)
    ka, kb = data['query'].map(seq_key), data['text'].map(seq_key)
    data['overlap'] = [overlap(a, b, seen_pos, seen_neg, proteins) for a, b in zip(ka, kb)]
    data['train_pos_test_neg'] = [frozenset((a, b)) in seen_pos and lab == 0 for a, b, lab in zip(ka, kb, data.label)]
    data['la'], data['lb'] = data['query'].str.len(), data['text'].str.len()
    data['tokens'] = data.la + data.lb + 3
    data['group'] = np.where(data.tokens <= MAX_TOKENS, 'short', 'long')
    data['bin'] = pd.cut(data.tokens, BINS, labels=BIN_LABELS).astype(str)
    data['has_longer_than_training'] = data[['la', 'lb']].max(1) > longest
    data['primary'] = data.overlap != 'pair_seen'
    if args.limit:                                          # smoke test only: a balanced random subset
        keep = data[data.primary].groupby('label', group_keys=False).apply(
            lambda g: g.sample(n=min(args.limit // 2, len(g)), random_state=SEED)).uid
        data = data[data.uid.isin(set(keep))].reset_index(drop=True)
    primary = data[data.primary]
    rng = np.random.default_rng(SEED)
    e2_ids = e2_sample(primary, args.e2_n, rng)
    data['e2'] = data.uid.isin(e2_ids)
    seqs = dict(zip(data.Uniprot_a, data['query'])); seqs.update(zip(data.Uniprot_b, data['text']))
    full = pd.read_csv(test, dtype=str, usecols=['Uniprot_a', 'Uniprot_b', 'label', 'query', 'text'])
    seqs_all = dict(zip(full.Uniprot_a, full['query'].str.strip().str.upper()))
    seqs_all.update(zip(full.Uniprot_b, full['text'].str.strip().str.upper()))
    partners = {}
    for a, b in zip(full.Uniprot_a[full.label == '1'], full.Uniprot_b[full.label == '1']):
        partners.setdefault(a, set()).add(b); partners.setdefault(b, set()).add(a)
    donors = sorted(seqs_all)
    jobs, e2_rows = [], []
    for row in primary.itertuples(index=False):
        jobs += [(f'{row.uid}|orig|AB', row.query, row.text), (f'{row.uid}|orig|BA', row.text, row.query)]
    for number, row in enumerate(data[data.e2].sort_values('uid').itertuples(index=False)):
        need = max(TARGETS) - row.tokens
        for pad in PADS:
            pa, pb = pad_streams(row, seqs_all, partners, donors, pad, need // 2, need - need // 2, number)
            for T in TARGETS:
                add = T - row.tokens
                qa, qb = row.query + pa[:add // 2], row.text + pb[:add - add // 2]
                check(len(qa) + len(qb) + 3 == T and qa.startswith(row.query) and qb.startswith(row.text),
                      f'{row.uid}: padding error')
                variant = f'{pad}{T}'
                jobs += [(f'{row.uid}|{variant}|AB', qa, qb), (f'{row.uid}|{variant}|BA', qb, qa)]
                e2_rows.append(dict(uid=row.uid, variant=variant, pad=pad, target=T, added=add))
    jobs = pd.DataFrame(jobs, columns=['pair_id', 'query', 'text'])
    check(jobs.pair_id.is_unique, 'Duplicate job ids')
    folder = out/'jobs'; folder.mkdir()
    for c, start in enumerate(range(0, len(jobs), args.chunk)):
        jobs.iloc[start:start + args.chunk].to_csv(folder/f'chunk_{c:04d}.csv', index=False)
    cols = ['uid', 'Uniprot_a', 'Uniprot_b', 'label', 'la', 'lb', 'tokens', 'group', 'bin', 'overlap',
            'has_longer_than_training', 'primary', 'e2']
    data[cols].to_csv(out/'pairs.csv', index=False)
    pd.DataFrame(e2_rows).to_csv(out/'e2_variants.csv', index=False)
    p = data[data.primary]
    summary = dict(created_utc=now(), test_file=str(test.resolve()), test_sha256=TEST_SHA256 if not args.skip_sha else digest(test),
        training_file=str(train.resolve()), training_sha256=TRAIN_SHA256 if not args.skip_sha else digest(train),
        bernett_revision=BERNETT_REV, longest_training_protein_aa=longest, limit=args.limit,
        pairs=int(len(data)), excluded_pair_seen=int((~data.primary).sum()),
        excluded_train_pos_test_neg=int(data.train_pos_test_neg.sum()),
        primary=dict(pairs=int(len(p)), positives=int((p.label == 1).sum()), negatives=int((p.label == 0).sum())),
        by_group=counts(p, ['group', 'label']),
        by_bin=counts(p, ['bin', 'label']),
        by_overlap=counts(p, ['group', 'overlap', 'label']),
        longer_than_training={str(k): int(v) for k, v in p.groupby('label').has_longer_than_training.sum().items()},
        e2=dict(pairs=int(data.e2.sum()), positives=int((data.e2 & (data.label == 1)).sum()),
                by_overlap=counts(data[data.e2], ['overlap', 'label'])),
        jobs=int(len(jobs)), chunks=int(np.ceil(len(jobs) / args.chunk)), chunk_size=args.chunk,
        max_job_tokens=int((jobs['query'].str.len() + jobs['text'].str.len() + 3).max()))
    write_json(out/'dataset.json', summary)
    print(json.dumps({k: summary[k] for k in ['pairs', 'excluded_pair_seen', 'primary', 'e2', 'jobs', 'chunks',
                                              'longest_training_protein_aa']}, indent=1))


# ------------------------------------------------------------------ work / merge
def work(args):
    runs = Path(args.runs)
    chunks = sorted((runs/'data'/'jobs').glob('chunk_*.csv'))
    check(chunks, f'No chunks under {runs}/data/jobs')
    for i, chunk in enumerate(chunks):
        if i % args.of != args.worker:
            continue
        target = runs/'scores'/chunk.name
        if Path(str(target) + '.meta.json').exists():
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists():                                   # interrupted: plm_score.py refuses to overwrite
            target.rename(target.with_suffix(f'.incomplete_{datetime.now():%H%M%S}'))
        cmd = [args.plm_python, str(ROOT/'plm_score.py'), 'score', '--jobs', str(chunk), '--base-model', args.base_model,
               '--checkpoint', args.checkpoint, '--checkpoint-name', args.checkpoint_name,
               '--checkpoint-revision', args.checkpoint_revision, '--allow-long', '--output', str(target)]
        if args.device:
            cmd += ['--device', args.device]
        print(f'{now()} worker {args.worker}/{args.of}: {chunk.name}', flush=True)
        with target.with_suffix('.log').open('w') as log:
            subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    print(f'{now()} worker {args.worker}/{args.of}: done', flush=True)


def merge(args):
    runs = Path(args.runs)
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    chunks = sorted((runs/'data'/'jobs').glob('chunk_*.csv'))
    frames, metas = [], []
    for chunk in chunks:
        target = runs/'scores'/chunk.name
        meta = Path(str(target) + '.meta.json')
        check(meta.exists(), f'{chunk.name} not scored')
        m = json.loads(meta.read_text())
        check(m.get('status') == 'complete' and m.get('input_sha256') == digest(chunk), f'{chunk.name}: meta mismatch')
        s = pd.read_csv(target)
        j = pd.read_csv(chunk, dtype=str)
        check(len(s) == len(j) and list(s.pair_id) == list(j.pair_id), f'{chunk.name}: jobs and scores differ')
        check((s.n_tokens.to_numpy() == (j['query'].str.len() + j['text'].str.len() + 3).to_numpy()).all(),
              f'{chunk.name}: token counts differ (truncation?)')
        frames.append(s[['pair_id', 'score']]); metas.append(m)
    s = pd.concat(frames, ignore_index=True)
    check(s.pair_id.is_unique and s.score.notna().all(), 'Duplicate or missing scores')
    parts = s.pair_id.str.split('|', expand=True)
    s['uid'], s['variant'], s['order'] = parts[0], parts[1], parts[2]
    w = s.pivot_table(index=['uid', 'variant'], columns='order', values='score', aggfunc='first').reset_index()
    check(w[['AB', 'BA']].notna().all().all(), 'A pair lacks one input order')
    w['score'] = (w.AB + w.BA) / 2
    w.rename(columns={'AB': 'score_ab', 'BA': 'score_ba'})[['uid', 'variant', 'score', 'score_ab', 'score_ba']].to_csv(out, index=False)
    ck = {m.get('checkpoint_sha256') for m in metas}
    write_json(Path(str(out) + '.json'), dict(created_utc=now(), chunks=len(chunks), jobs=int(len(s)),
        pairs_variants=int(len(w)), checkpoint_sha256=sorted(x for x in ck if x), output_sha256=digest(out)))
    print(f'{len(w)} pair-variants -> {out}')


# ------------------------------------------------------------------ statistics
def auroc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if not len(pos) or not len(neg):
        return np.nan
    r = rankdata(np.concatenate([pos, neg]))
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def aupr(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    scores = np.concatenate([pos, neg]); labels = np.r_[np.ones(len(pos)), np.zeros(len(neg))]
    order = np.argsort(-scores, kind='stable'); scores, labels = scores[order], labels[order]
    last = np.r_[np.flatnonzero(np.diff(scores)), len(scores) - 1]
    tp = np.cumsum(labels)[last]; fp = (last + 1) - tp
    return float(np.sum(np.diff(np.r_[0, tp / len(pos)]) * tp / (tp + fp)))


def ci(values):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    return [float(x) for x in np.percentile(v, [2.5, 97.5])] if len(v) else [np.nan, np.nan]


def delta_groups(t, rng, n_boot=N_BOOT):
    """AUROC(short) - AUROC(long); bootstrap within each group x label."""
    parts = {(g, l): t.score[(t.group == g) & (t.label == l)].to_numpy() for g in ('short', 'long') for l in (1, 0)}
    a_s, a_l = auroc(parts['short', 1], parts['short', 0]), auroc(parts['long', 1], parts['long', 0])
    boots = []
    for _ in range(n_boot):
        r = {k: v[rng.integers(len(v), size=len(v))] if len(v) else v for k, v in parts.items()}
        boots.append(auroc(r['short', 1], r['short', 0]) - auroc(r['long', 1], r['long', 0]))
    return dict(auroc_short=a_s, auroc_long=a_l, delta=a_s - a_l, ci95=ci(boots),
                n={f'{g}_{"pos" if l else "neg"}': int(len(v)) for (g, l), v in parts.items()})


def auroc_ci(pos, neg, rng, n_boot=N_BOOT):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    boots = [auroc(pos[rng.integers(len(pos), size=len(pos))], neg[rng.integers(len(neg), size=len(neg))])
             for _ in range(n_boot)] if len(pos) and len(neg) else []
    return dict(auroc=auroc(pos, neg), ci95=ci(boots), n_pos=int(len(pos)), n_neg=int(len(neg)))


def protein_splits(t, n_splits=N_SPLITS):
    proteins = np.array(sorted(set(t.Uniprot_a) | set(t.Uniprot_b)))
    out = []
    for seed in range(1, n_splits + 1):
        half = dict(zip(proteins, np.random.default_rng(seed).integers(2, size=len(proteins))))
        ha, hb = t.Uniprot_a.map(half), t.Uniprot_b.map(half)
        for h in (0, 1):
            sub = t[(ha == h) & (hb == h)]
            s = lambda g, l: sub.score[(sub.group == g) & (sub.label == l)]
            out.append(dict(seed=seed, half=h, pairs=int(len(sub)),
                            delta=auroc(s('short', 1), s('short', 0)) - auroc(s('long', 1), s('long', 0))))
    return pd.DataFrame(out)


def e1_decision(overall, strata):
    d, lo, hi = overall['delta'], *overall['ci95']
    signs = [np.sign(strata[k]['delta']) for k in STRATA if np.isfinite(strata[k]['delta'])]
    same = len(signs) == len(STRATA) and len(set(signs)) == 1 and signs[0] == np.sign(d)
    if d >= MARGIN and lo > 0 and same:
        return 'long_clearly_worse'
    if lo > 0:
        return 'long_slightly_worse'
    if -MARGIN <= lo and hi <= MARGIN:
        return 'no_material_difference'
    return 'uncertain'


def e2_table(t, rng, n_boot=N_BOOT):
    """AUROC original vs padded on the same E2 pairs; paired bootstrap over positives and negatives."""
    base = t[t.variant == 'orig'].set_index('uid')
    rows = []
    for variant in sorted(set(t.variant) - {'orig'}):
        pad = t[t.variant == variant].set_index('uid')
        ids = base.index.intersection(pad.index)
        lab = base.loc[ids, 'label'].to_numpy()
        so, sp = base.loc[ids, 'score'].to_numpy(), pad.loc[ids, 'score'].to_numpy()
        P, N = np.flatnonzero(lab == 1), np.flatnonzero(lab == 0)
        a0, a1 = auroc(so[P], so[N]), auroc(sp[P], sp[N])
        boots = []
        for _ in range(n_boot):
            p, n = P[rng.integers(len(P), size=len(P))], N[rng.integers(len(N), size=len(N))]
            boots.append(auroc(so[p], so[n]) - auroc(sp[p], sp[n]))
        m = re.fullmatch(r'([a-z]+)(\d+)', variant)
        check(m is not None, f'Unexpected variant {variant}')
        rows.append(dict(variant=variant, pad=m.group(1), target=int(m.group(2)),
                         n_pos=int(len(P)), n_neg=int(len(N)), auroc_original=a0, auroc_padded=a1, delta=a0 - a1,
                         ci95_low=ci(boots)[0], ci95_high=ci(boots)[1],
                         median_pos_original=float(np.median(so[P])), median_pos_padded=float(np.median(sp[P]))))
    return pd.DataFrame(rows)


def e2_decision(table, target=2400):
    rows = table[table.target == target]
    if len(rows) != len(PADS):
        return 'uncertain'
    if ((rows.delta >= MARGIN) & (rows.ci95_low > 0)).all():
        return 'length_itself_lowers_discrimination'
    if ((rows.ci95_low >= -MARGIN) & (rows.ci95_high <= MARGIN)).all():
        return 'length_itself_small_effect'
    return 'uncertain'


CONCLUSION = {
    ('worse', 'length_itself_lowers_discrimination'): 'V12 的长序列问题来自模型对长度的处理;胞外段打分或裁剪值得做',
    ('worse', 'length_itself_small_effect'): '长配对本身更难(序列内容或标签),不是长度造成的;裁剪不一定有用',
    ('none', None): '长度不是 V12 的主要问题;小鼠面板上的长度差异要另找原因',
}


def conclusion(e1, e2):
    if e1 in ('long_clearly_worse', 'long_slightly_worse') and e2 in ('length_itself_lowers_discrimination',
                                                                     'length_itself_small_effect'):
        return CONCLUSION['worse', e2]
    if e1 == 'no_material_difference':
        return CONCLUSION['none', None]
    return '其他组合:如实报告'


def lr_replication(path, rng):
    """H-A humanV12 scores (plm_lr_benchmark): single-gene, non-homophilic positives vs random negatives,
    pairs seen in V12 training removed; same token rule. Descriptive."""
    t = pd.read_csv(path)
    t = t[t.negative.isin(['positive', 'random']) & ~t.training.isin(['seen_positive', 'seen_negative'])].copy()
    t['tokens'] = pd.to_numeric(t.ligand_lengths) + pd.to_numeric(t.receptor_lengths) + 3
    t['group'] = np.where(t.tokens <= MAX_TOKENS, 'short', 'long')
    t['label'] = (t.negative == 'positive').astype(int)
    res = dict(overall=delta_groups(t, rng))
    for cat in ['proteins_seen', 'protein_unseen']:
        res[cat] = delta_groups(t[t.training == cat], rng)
    return res


def plot(bins, e2, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
    ax = axes[0]
    x = np.arange(len(bins))
    ax.errorbar(x, bins.auroc, yerr=[bins.auroc - bins.ci95_low, bins.ci95_high - bins.auroc], fmt='o-', color='#1f5fa8')
    ax.axvline(2.5, color='0.6', lw=0.8, ls='--'); ax.text(2.55, ax.get_ylim()[0], ' 1603 tokens', color='0.4', fontsize=7, va='bottom')
    ax.set_xticks(x); ax.set_xticklabels(bins.bin, fontsize=7); ax.set_xlabel('Tokens (both sequences)'); ax.set_ylabel('AUROC')
    ax.set_title('E1: AUROC by length, primary set', loc='left', fontsize=9)
    ax = axes[1]
    for pad, color in zip(PADS, ['#1f5fa8', '#c2571a']):
        r = e2[e2['pad'] == pad].sort_values('target')
        xs = [E2_MAX_TOKENS] + list(r.target); ys = [r.auroc_original.iloc[0]] + list(r.auroc_padded)
        ax.plot(xs, ys, 'o-', color=color, label=f'padding: {pad}')
    ax.set_xlabel('Tokens after padding (original pairs <= 1000)'); ax.set_ylabel('AUROC'); ax.legend(frameon=False, fontsize=7)
    ax.set_title('E2: same pairs, padded', loc='left', fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=200); plt.close(fig)


def evaluate(args):
    out = fresh(args.output_dir)
    pairs = pd.read_csv(Path(args.runs)/'data'/'pairs.csv', dtype={'uid': str})
    scores = pd.read_csv(args.scores, dtype={'uid': str})
    t = scores.merge(pairs, on='uid', how='left')
    check(t.label.notna().all(), 'Scores for unknown pairs')
    prim = t[(t.variant == 'orig') & t.primary]
    check(len(prim) == int(pairs.primary.sum()), 'Not every primary pair has a score')
    rng = np.random.default_rng(SEED)
    overall = delta_groups(prim, rng, args.n_boot)
    strata = {k: delta_groups(prim[prim.overlap == k], rng, args.n_boot) for k in STRATA}
    splits = protein_splits(prim)
    finite = splits.delta[np.isfinite(splits.delta)]
    split_same = float((np.sign(finite) == np.sign(overall['delta'])).mean()) if len(finite) else np.nan
    bins = []
    for b in BIN_LABELS:
        s = prim[prim.bin == b]
        r = auroc_ci(s.score[s.label == 1], s.score[s.label == 0], rng, args.n_boot)
        bins.append(dict(bin=b, auroc=r['auroc'], ci95_low=r['ci95'][0], ci95_high=r['ci95'][1], n_pos=r['n_pos'], n_neg=r['n_neg']))
    bins = pd.DataFrame(bins)
    longp = prim[prim.group == 'long']
    beyond = {str(k): auroc_ci(g.score[g.label == 1], g.score[g.label == 0], rng, args.n_boot)
              for k, g in longp.groupby('has_longer_than_training')}
    descr = {g: dict(aupr=aupr(prim.score[(prim.group == g) & (prim.label == 1)], prim.score[(prim.group == g) & (prim.label == 0)]),
                     median_pos=float(prim.score[(prim.group == g) & (prim.label == 1)].median()),
                     median_neg=float(prim.score[(prim.group == g) & (prim.label == 0)].median())) for g in ('short', 'long')}
    e1 = e1_decision(overall, strata)
    e2t = e2_table(t[t.e2], rng, args.n_boot)
    e2 = e2_decision(e2t)
    lr = lr_replication(args.lr_scores, rng) if args.lr_scores else None
    metrics = dict(created_utc=now(), n_boot=args.n_boot, margin=MARGIN, E1=dict(overall=overall, strata=strata,
                   protein_splits=dict(n=int(len(finite)), same_sign_fraction=split_same,
                                       delta_min=float(finite.min()) if len(finite) else None,
                                       delta_max=float(finite.max()) if len(finite) else None),
                   longer_than_training_within_long=beyond, descriptive=descr),
                   E2=e2t.to_dict('records'), LR_replication=lr)
    verdict = dict(E1=e1, E2=e2, conclusion=conclusion(e1, e2),
                   E1_delta=overall['delta'], E1_ci95=overall['ci95'],
                   E2_at_2400={r['pad']: dict(delta=r['delta'], ci95=[r['ci95_low'], r['ci95_high']])
                               for r in e2t[e2t.target == 2400].to_dict('records')})
    write_json(out/'metrics.json', metrics); write_json(out/'verdict.json', verdict)
    bins.to_csv(out/'e1_bins.csv', index=False); splits.to_csv(out/'e1_protein_splits.csv', index=False)
    e2t.to_csv(out/'e2_padding.csv', index=False)
    if len(e2t):
        plot(bins, e2t, out/'length_effect.png')
    print(json.dumps(verdict, indent=2, ensure_ascii=False, default=float))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fetch'); f.add_argument('--output', required=True)
    q = sub.add_parser('prepare')
    q.add_argument('--test-file', required=True); q.add_argument('--training-root', required=True)
    q.add_argument('--output-dir', required=True)
    q.add_argument('--limit', type=int, default=0, help='smoke test only: balanced random subset of the primary set')
    q.add_argument('--e2-n', type=int, default=E2_N, help='positives and negatives each (plan: 2000)')
    q.add_argument('--chunk', type=int, default=CHUNK)
    q.add_argument('--skip-sha', action='store_true', help='unit tests only')
    w = sub.add_parser('work')
    w.add_argument('--runs', required=True); w.add_argument('--worker', type=int, default=0); w.add_argument('--of', type=int, default=1)
    w.add_argument('--plm-python', required=True, help='python of the PLM-interact environment')
    w.add_argument('--checkpoint', required=True); w.add_argument('--base-model', required=True)
    w.add_argument('--checkpoint-name', default=CHECKPOINT[0]); w.add_argument('--checkpoint-revision', default=CHECKPOINT[1])
    w.add_argument('--device')
    m = sub.add_parser('merge'); m.add_argument('--runs', required=True); m.add_argument('--output', required=True)
    e = sub.add_parser('evaluate')
    e.add_argument('--runs', required=True); e.add_argument('--scores', required=True); e.add_argument('--output-dir', required=True)
    e.add_argument('--lr-scores', help='reports/plm_lr_benchmark/results/pair_scores/H-A__humanV12.csv.gz')
    e.add_argument('--n-boot', type=int, default=N_BOOT)
    return p


if __name__ == '__main__':
    a = parser().parse_args()
    dict(fetch=fetch, prepare=prepare, work=work, merge=merge, evaluate=evaluate)[a.command](a)
