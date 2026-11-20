"""Experiment 3 (reports/affinity_01/plan.md): can BindPred rank measured binding affinities
of endogenous peptide-receptor pairs well enough to serve as S in the OT cost?

  prepare   GtoPdb + UniProt + BindPred training table -> pairs.csv, excluded.csv, plm_input.csv
  embed     unique sequences -> embeddings.npz (frozen ESM-2 650M, BindPred recipe)
  predict   embeddings + CatBoost model -> bindpred.csv (log10 Kd)
  evaluate  pairs + scores -> summary/{decisions.json, groups.csv, pairs_scored.csv, figure}

Environments: `prepare`/`evaluate` need pandas+scipy+biopython (commot-research);
`embed` needs torch+transformers (the PLM env); `predict` needs catboost.

BindPred embedding recipe, verified against its published training table (max abs diff 3e-6):
mean over residue token hidden states of the last layer, sequences truncated to the first
1022 residues, ligand features first then receptor features (2560 dims). The model is
order-sensitive: ligand-first PCC 0.926 vs receptor-first 0.696 on its own training rows.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from research import check

MAX_RESIDUES = 1022                  # ESM-2 position limit used by the BindPred authors
MIN_GROUP = 3                        # plan section 3: groups need >= 3 members
MIN_GROUPS = 10                      # plan section 4: fewer usable groups -> undetermined
RHO_BAR_MIN = 0.30                   # plan section 4, A2
N_BOOT, BOOT_SEED = 1000, 1          # plan section 3
IDENTITY = 0.90                      # plan section 2, leakage containment threshold
AA = set('ACDEFGHIKLMNPQRSTVWY')
SCORERS = ('bindpred', 'v12', 'length')


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as s:
        for block in iter(lambda: s.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, default=lambda x: x.item() if hasattr(x, 'item') else str(x)))


def read_gtopdb(path, **kw):
    """GtoPdb CSVs carry a version comment on line 1."""
    first = Path(path).read_text(encoding='utf-8', errors='replace').split('\n', 1)[0]
    return pd.read_csv(path, skiprows=1, low_memory=False, **kw), first.strip().strip('"')


# ------------------------------------------------------------------ prepare
def candidate_pairs(gtopdb_dir, uniprot):
    """Plan section 2 inclusion rules. Returns (pairs, counts)."""
    gd = Path(gtopdb_dir)
    it, version = read_gtopdb(gd/'gtopdb_interactions.csv')
    pep, version_pep = read_gtopdb(gd/'gtopdb_peptides.csv')
    check(version == version_pep, 'GtoPdb files come from different releases')
    up = pd.read_csv(uniprot, sep='\t')
    seqs = dict(zip(up.accession, up.sequence))
    counts = dict(interactions=len(it))
    d = it[it['Endogenous'].astype(str).str.lower() == 'true']
    counts['endogenous'] = len(d)
    d = d[d['Ligand Type'] == 'Peptide']
    counts['peptide'] = len(d)
    d = d[d['Affinity Units'].isin(['pKd', 'pKi'])]
    counts['binding_units'] = len(d)
    d = d[~d['Target UniProt ID'].astype(str).str.contains(r'\|', regex=True)]
    counts['single_gene_target'] = len(d)
    d = d.copy()
    d['value'] = d['Affinity Median'].where(d['Affinity Median'].notna(),
                                            (d['Affinity High'] + d['Affinity Low']) / 2)
    d = d[d['value'].notna()]
    counts['with_value'] = len(d)
    seq_col = 'Single letter amino acid sequence'
    pep = pep[['Ligand id', seq_col]].rename(columns={'Ligand id': 'lig_id'})
    pep['lig_id'] = pep['lig_id'].astype(str)
    d['lig_id'] = d['Ligand ID'].astype(str)
    d = d.merge(pep, on='lig_id', how='left')
    raw = d[seq_col].astype('string').fillna('').str.upper().astype(str)   # pandas keeps NA through astype(str)
    d['pyroglutamate'] = raw.str.startswith('X')
    d['lig_seq'] = raw.str.replace('X', 'Q', regex=False)
    ok = d['lig_seq'].map(lambda s: len(s) >= 5 and set(s) <= AA)
    counts['ligand_sequence_ok'] = int(ok.sum())
    counts['pyroglutamate_renamed'] = int(d.loc[ok, 'pyroglutamate'].sum())
    d = d[ok]
    d['rec_acc'] = d['Target UniProt ID'].astype(str)
    d['rec_seq'] = d['rec_acc'].map(seqs)
    check(d['rec_seq'].notna().all(), 'A target accession is missing from the UniProt file')
    g = d.groupby(['Target Species', 'rec_acc', 'lig_id'], as_index=False).agg(
        target=('Target', 'first'), ligand=('Ligand', 'first'), rec_seq=('rec_seq', 'first'),
        lig_seq=('lig_seq', 'first'), value=('value', 'median'), n_records=('value', 'size'),
        units=('Affinity Units', lambda x: ','.join(sorted(set(x)))))
    g['species'] = g['Target Species']
    g['pair_id'] = g.species.str.lower() + '__' + g.rec_acc + '__' + g.lig_id
    g['rec_len'] = g.rec_seq.str.len()
    g['lig_len'] = g.lig_seq.str.len()
    counts['pairs'] = len(g)
    counts['gtopdb_version'] = version
    return g.drop(columns='Target Species').reset_index(drop=True), counts


def leakage(pairs, training_csv, identity=IDENTITY):
    """Plan section 2: drop pairs whose ligand and receptor both match the two sequences of one
    training row (>= `identity` of the shorter sequence aligned identically)."""
    from collections import defaultdict
    from Bio import Align
    from Bio.Align import substitution_matrices
    tr = pd.read_csv(training_csv, usecols=['Receptor Sequence', 'Ligand Sequence']).dropna()
    seqs = sorted(set(tr['Receptor Sequence']) | set(tr['Ligand Sequence']))
    at = {s: i for i, s in enumerate(seqs)}
    rows = {(at[r], at[l]) for r, l in zip(tr['Receptor Sequence'], tr['Ligand Sequence'])}
    rows |= {(b, a) for a, b in rows}
    k = 4
    index = defaultdict(set)
    n_kmers = []
    for j, s in enumerate(seqs):
        kmers = {s[i:i+k] for i in range(len(s) - k + 1)}
        n_kmers.append(len(kmers))
        for kmer in kmers:
            index[kmer].add(j)
    aligner = Align.PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                                    open_gap_score=-10, extend_gap_score=-0.5)
    def matches(q):
        kmers = {q[i:i+k] for i in range(len(q) - k + 1)}
        counted = defaultdict(int)
        for kmer in kmers:
            for j in index.get(kmer, ()):
                counted[j] += 1
        out = set()
        for j, c in counted.items():
            t = seqs[j]
            shorter = min(len(q), len(t))
            if c < max(1, 0.5 * min(len(kmers), n_kmers[j])):
                continue      # cheap filter before the alignment: half the smaller k-mer set must be shared
            a = aligner.align(q, t)[0]
            same = sum(1 for (qs, qe), (ts, te) in zip(*a.aligned)
                       for x, y in zip(q[qs:qe], t[ts:te]) if x == y)
            if same >= identity * shorter:
                out.add(j)
        return out
    hits = {q: matches(q) for q in sorted(set(pairs.rec_seq) | set(pairs.lig_seq))}
    p = pairs.copy()
    p['lig_seen'] = [bool(hits[s]) for s in p.lig_seq]
    p['rec_seen'] = [bool(hits[s]) for s in p.rec_seq]
    p['pair_seen'] = [any((a, b) in rows for a in hits[l] for b in hits[r])
                      for l, r in zip(p.lig_seq, p.rec_seq)]
    return p


def prepare(args):
    out = Path(args.output_dir)
    check(not out.exists(), f'{out} exists')
    pairs, counts = candidate_pairs(args.gtopdb_dir, args.uniprot)
    pairs = leakage(pairs, args.training_csv)
    counts['excluded_pair_seen'] = int(pairs.pair_seen.sum())
    kept = pairs[~pairs.pair_seen].reset_index(drop=True)
    if args.limit:
        groups = kept.groupby(['species', 'rec_acc']).pair_id.count().sort_values(ascending=False)
        keep = {g for g, _ in list(groups.items())[:args.limit]}
        kept = kept[[(s, r) in keep for s, r in zip(kept.species, kept.rec_acc)]].reset_index(drop=True)
    counts['primary'] = len(kept)
    counts['ligand_seen_only'] = int((kept.lig_seen & ~kept.rec_seen).sum())
    counts['receptor_seen_only'] = int((kept.rec_seen & ~kept.lig_seen).sum())
    counts['receptors_truncated'] = int((kept.rec_len > MAX_RESIDUES).sum())
    counts['ligands_truncated'] = int((kept.lig_len > MAX_RESIDUES).sum())
    out.mkdir(parents=True)
    kept.to_csv(out/'pairs.csv', index=False)
    pairs[pairs.pair_seen].to_csv(out/'excluded.csv', index=False)
    kept.rename(columns={'lig_seq': 'query', 'rec_seq': 'text'})[['pair_id', 'query', 'text']] \
        .to_csv(out/'plm_input.csv', index=False)
    write_json(out/'dataset.json', dict(status='complete', created_utc=now(), limit=args.limit, counts=counts,
               inputs={f: digest(Path(args.gtopdb_dir)/f) for f in ['gtopdb_interactions.csv', 'gtopdb_peptides.csv']}
               | {'uniprot': digest(args.uniprot), 'bindpred_training': digest(args.training_csv)},
               outputs={f: digest(out/f) for f in ['pairs.csv', 'excluded.csv', 'plm_input.csv']}))
    print(json.dumps(counts, indent=2))


# ------------------------------------------------------------------ embed
def embed(args):
    """Frozen ESM-2 650M, mean over residue tokens, first MAX_RESIDUES residues (BindPred recipe)."""
    import torch
    from transformers import AutoTokenizer, EsmModel
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    pairs = pd.read_csv(args.pairs)
    seqs = sorted(set(pairs.lig_seq) | set(pairs.rec_seq))
    torch.set_num_threads(args.threads)
    tok = AutoTokenizer.from_pretrained(args.esm2)
    model = EsmModel.from_pretrained(args.esm2, add_pooling_layer=False).eval()
    device = torch.device(args.device or ('cuda:0' if torch.cuda.is_available() else 'cpu'))
    model.to(device)
    vecs = np.zeros((len(seqs), model.config.hidden_size), dtype=np.float64)
    for i, s in enumerate(seqs):
        ids = tok(s[:MAX_RESIDUES], return_tensors='pt', add_special_tokens=True)
        n = len(s[:MAX_RESIDUES])
        check(int(ids['input_ids'].shape[1]) == n + 2, f'Unexpected token count for sequence {i}')
        with torch.no_grad():
            h = model(**{k: v.to(device) for k, v in ids.items()}).last_hidden_state[0]
        vecs[i] = h[1:n+1].double().mean(0).cpu().numpy()
        if (i + 1) % 25 == 0:
            print(f'{i+1}/{len(seqs)}', flush=True)
    np.savez_compressed(out, sequences=np.asarray(seqs, object), vectors=vecs,
                        meta=json.dumps(dict(esm2=str(Path(args.esm2).resolve()), device=str(device),
                                             torch=torch.__version__, created_utc=now(),
                                             hidden_size=model.config.hidden_size, max_residues=MAX_RESIDUES)))
    print(f'Embedded {len(seqs)} sequences -> {out}')


# ------------------------------------------------------------------ predict
def predict(args):
    from catboost import CatBoostRegressor
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    pairs = pd.read_csv(args.pairs)
    z = np.load(args.embeddings, allow_pickle=True)
    vec = dict(zip(z['sequences'].tolist(), z['vectors']))
    missing = {s for s in set(pairs.lig_seq) | set(pairs.rec_seq) if s not in vec}
    check(not missing, f'{len(missing)} sequences are missing from the embedding file')
    X = np.vstack([np.r_[vec[l], vec[r]] for l, r in zip(pairs.lig_seq, pairs.rec_seq)])   # ligand first
    model = CatBoostRegressor()
    model.load_model(args.model, format='cbm')
    log10_kd = model.predict(X)
    pd.DataFrame(dict(pair_id=pairs.pair_id, log10_kd_pred=log10_kd, pkd_hat=-log10_kd)).to_csv(out, index=False)
    write_json(str(out) + '.meta.json', dict(status='complete', created_utc=now(), n_pairs=len(pairs),
               model=str(Path(args.model).resolve()), model_sha256=digest(args.model),
               embeddings_sha256=digest(args.embeddings), pairs_sha256=digest(args.pairs),
               trees=int(model.tree_count_), features=int(X.shape[1]), output_sha256=None))
    print(f'Predicted {len(pairs)} pairs -> {out}')


# ------------------------------------------------------------------ evaluate
def scored_table(args):
    """pairs.csv + predictions -> one row per pair with every scorer's value (plan section 1)."""
    pairs = pd.read_csv(args.pairs)
    bp = pd.read_csv(args.bindpred)
    t = pairs.merge(bp, on='pair_id', how='left', validate='one_to_one')
    check(t.pkd_hat.notna().all(), 'A pair has no BindPred prediction')
    v12 = pd.read_csv(args.v12)
    check(v12.pair_id.is_unique, 'Duplicate pair_id in the V12 scores')
    t = t.merge(v12[['pair_id', 'score']].rename(columns={'score': 'v12'}), on='pair_id', how='left',
                validate='one_to_one')
    check(t.v12.notna().all(), 'A pair has no V12 score')
    t['bindpred'] = t.pkd_hat
    t['length'] = t.lig_len.astype(float)
    t['log10_kd_exp'] = -t.value
    return t


def group_rho(table, key, member, scorers=SCORERS, min_group=MIN_GROUP):
    """Within-group Spearman for each scorer; one row per usable group (plan section 3)."""
    from scipy.stats import spearmanr
    rows = []
    for k, g in table.groupby(key):
        g = g.drop_duplicates(member)
        if len(g) < min_group or g.value.nunique() < 2:
            continue                                        # ties only: Spearman undefined
        r = dict(group='__'.join(map(str, k)) if isinstance(k, tuple) else str(k), n=len(g),
                 value_range=float(g.value.max() - g.value.min()),
                 informative_pairs=int(sum(abs(a - b) >= 0.5 for i, a in enumerate(g.value)
                                           for b in g.value.to_numpy()[i+1:])))
        for s in scorers:
            r[s] = float(spearmanr(g[s], g.value)[0]) if g[s].nunique() > 1 else np.nan
        rows.append(r)
    return pd.DataFrame(rows)


def rho_bars(groups, scorers=SCORERS):
    """Mean within-group Spearman. The length baseline takes the favourable sign (plan section 4)."""
    out = {}
    for s in scorers:
        v = groups[s].to_numpy(float)
        m = float(np.nanmean(v)) if np.isfinite(v).any() else float('nan')
        out[s] = abs(m) if s == 'length' else m
    return out


def bootstrap(groups, scorers=SCORERS, n_boot=N_BOOT, seed=BOOT_SEED):
    """Resample groups (plan section 3); returns per-scorer draws and paired differences."""
    rng = np.random.default_rng(seed)
    draws = {s: [] for s in scorers}
    diffs = {s: [] for s in scorers if s != 'bindpred'}
    idx = np.arange(len(groups))
    for _ in range(n_boot):
        take = groups.iloc[rng.choice(idx, len(idx), replace=True)]
        bars = rho_bars(take, scorers)
        for s in scorers:
            draws[s].append(bars[s])
        for s in diffs:
            diffs[s].append(bars['bindpred'] - bars[s])
    ci = lambda v: [float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))]
    return {s: ci(v) for s, v in draws.items()}, {s: ci(v) for s, v in diffs.items()}


def decide(groups):
    """Plan section 4."""
    if len(groups) < MIN_GROUPS:
        return dict(undetermined=True, note=f'only {len(groups)} usable groups (< {MIN_GROUPS})',
                    n_groups=int(len(groups)))
    bars = rho_bars(groups)
    ci, diff = bootstrap(groups)
    a1 = bool(ci['bindpred'][0] > 0)
    a2 = bool(bars['bindpred'] >= RHO_BAR_MIN)
    a3 = bool(all(diff[s][0] > 0 for s in diff))
    verdict = ('可作为 S 进入 OT' if a1 and a2 and a3 else
               '在本测评集上看不出排序能力' if not a1 else
               '有排序能力但效应量不足' if not a2 else '效应量足够但不优于基线')
    return dict(undetermined=False, n_groups=int(len(groups)), n_pairs_in_groups=int(groups.n.sum()),
                rho_bar=bars, ci95=ci, paired_diff_ci95=diff,
                A1_above_zero=a1, A2_effect_size=a2, A3_above_baselines=a3,
                rho_bar_min=RHO_BAR_MIN, verdict=verdict)


def strata(table, key, member):
    """Plan section 5.2, descriptive."""
    out = {}
    bins = [(0, 20), (21, 50), (51, 10 ** 6)]
    for name, sub in [('species_' + s, table[table.species == s]) for s in sorted(table.species.unique())] + \
                     [('ligand_seen_only', table[table.lig_seen & ~table.rec_seen]),
                      ('receptor_seen_only', table[table.rec_seen & ~table.lig_seen]),
                      ('neither_seen', table[~table.lig_seen & ~table.rec_seen]),
                      ('pKd_only', table[table.units == 'pKd']), ('pKi_only', table[table.units == 'pKi'])] + \
                     [(f'ligand_len_{a}_{b}', table[table.lig_len.between(a, b)]) for a, b in bins]:
        g = group_rho(sub, key, member)
        out[name] = dict(n_pairs=int(len(sub)), n_groups=int(len(g)),
                         rho_bar=rho_bars(g) if len(g) else None)
    return out


def figure(table, groups, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    ax = axes[0]
    for s, m in zip(sorted(table.species.unique()), 'osD^'):
        t = table[table.species == s]
        ax.scatter(t.log10_kd_exp, t.log10_kd_pred, s=22, marker=m, alpha=0.75, label=f'{s} (n={len(t)})')
    lo = float(min(table.log10_kd_exp.min(), table.log10_kd_pred.min())) - 0.5
    hi = float(max(table.log10_kd_exp.max(), table.log10_kd_pred.max())) + 0.5
    ax.plot([lo, hi], [lo, hi], color='0.3', linestyle='--', linewidth=1)
    ax.set(xlabel='experimental log10 Kd', ylabel='BindPred log10 Kd', xlim=(lo, hi), ylim=(lo, hi))
    ax.legend(frameon=False, fontsize=8, loc='upper left')
    ax = axes[1]
    data = [groups[s].dropna().to_numpy(float) for s in SCORERS]
    ax.axhline(0, color='0.6', linewidth=1)
    ax.boxplot(data, widths=0.55, showfliers=False)
    ax.set_xticks(range(1, len(SCORERS) + 1), list(SCORERS))
    rng = np.random.default_rng(0)
    for i, v in enumerate(data, start=1):
        ax.scatter(i + rng.uniform(-0.14, 0.14, len(v)), v, s=16, alpha=0.6)
    ax.set(ylabel='within-receptor Spearman', ylim=(-1.05, 1.05))
    fig.tight_layout()
    fig.savefig(path, dpi=200)


def evaluate(args):
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    table = scored_table(args)
    rec = group_rho(table, ['species', 'rec_acc'], 'lig_id')
    lig = group_rho(table, ['species', 'lig_id'], 'rec_acc')
    from scipy.stats import pearsonr, spearmanr
    err = table.log10_kd_pred - table.log10_kd_exp
    pooled = dict(n=int(len(table)), pearson=float(pearsonr(table.log10_kd_exp, table.log10_kd_pred)[0]),
                  spearman=float(spearmanr(table.log10_kd_exp, table.log10_kd_pred)[0]),
                  rmse=float(np.sqrt((err ** 2).mean())), bias=float(err.mean()),
                  v12_pooled_spearman=float(spearmanr(table.value, table.v12)[0]))
    result = dict(created_utc=now(), n_boot=N_BOOT, boot_seed=BOOT_SEED,
                  primary_receptor_groups=decide(rec),
                  secondary_ligand_groups=decide(lig) if len(lig) else dict(undetermined=True, note='no groups'),
                  pooled=pooled, strata=strata(table, ['species', 'rec_acc'], 'lig_id'),
                  group_summary=dict(median_value_range=float(rec.value_range.median()) if len(rec) else None,
                                     informative_pairs=int(rec.informative_pairs.sum()) if len(rec) else 0),
                  inputs={k: digest(getattr(args, k)) for k in ['pairs', 'bindpred', 'v12']})
    out.mkdir(parents=True)
    table.drop(columns=['rec_seq', 'lig_seq']).to_csv(out/'pairs_scored.csv', index=False)
    rec.to_csv(out/'groups_receptor.csv', index=False)
    lig.to_csv(out/'groups_ligand.csv', index=False)
    figure(table, rec, out/'affinity01_figure.png')
    write_json(out/'decisions.json', result)
    print(json.dumps({k: result[k] for k in ['primary_receptor_groups', 'pooled']}, indent=2, default=str))


def parser():
    pa = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = pa.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    p.add_argument('--gtopdb-dir', required=True, help='folder with gtopdb_interactions.csv and gtopdb_peptides.csv')
    p.add_argument('--uniprot', required=True); p.add_argument('--training-csv', required=True,
                   help='BindPred embeddings_only.csv (sequence columns are read, embeddings are not)')
    p.add_argument('--output-dir', required=True); p.add_argument('--limit', type=int, default=0,
                   help='smoke: keep only the N largest receptor groups')
    e = sub.add_parser('embed')
    e.add_argument('--pairs', required=True); e.add_argument('--esm2', required=True)
    e.add_argument('--output', required=True); e.add_argument('--device')
    e.add_argument('--threads', type=int, default=8)
    q = sub.add_parser('predict')
    q.add_argument('--pairs', required=True); q.add_argument('--embeddings', required=True)
    q.add_argument('--model', required=True); q.add_argument('--output', required=True)
    v = sub.add_parser('evaluate')
    for k in ('pairs', 'bindpred', 'v12', 'output'):
        v.add_argument('--' + k, required=True)
    return pa


if __name__ == '__main__':
    a = parser().parse_args()
    dict(prepare=prepare, embed=embed, predict=predict, evaluate=evaluate)[a.command](a)
