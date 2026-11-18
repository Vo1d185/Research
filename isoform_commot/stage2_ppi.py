"""PLM-interact scores for the stage-2 LR panel (step 1 of reports/stage2_ppi_01/plan.md). No OT here.

Environment: commot-research (numpy, pandas). Scoring itself: plm_score.py in the plminteract env.

  genes    stage-2 inputs (lr.csv, entities.csv) -> genes.csv (every subunit symbol) for
           `proteins.py map --genes genes.csv --species mouse ...`
  jobs     entity pairs -> unique subunit pairs -> plm_input.csv (both input orders, full length);
           optionally also random UNANNOTATED panel combinations (--cross-ratio), used only to estimate how
           well the scores separate panel pairs from unannotated combinations (AUROC); not negatives
  import   plm_score.py output -> ppi_scores.csv: one score per entity pair
             subunit pair score = mean of the two input orders
             entity pair score  = max over its subunit pairs (same rule as plm_lr_benchmark H-X);
                                  pairs of one protein with itself are skipped;
                                  defined only when every subunit is mapped (no imputation)
  designs  ppi_scores.csv + one slice's lr.csv -> designs_<slice>.csv with the score vectors that
           step 2 may use: raw S, S'(ligand) = S / sum of S over the receptors of that ligand,
           uniform 1/k(ligand), S'(receptor); plus descriptive statistics. Still no OT.

Every path is a command-line argument; existing outputs are never overwritten.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

OFFICIAL_MAX_TOKENS = 1603
EXPECTED_CHECKPOINT_SHA256 = '92565937d65a9b91a4e897c99653c972b99cfe37067eb42d9e783af688350547'   # humanV12
ORDERS = ('AB', 'BA')


def check(ok, msg):
    if not ok:
        raise ValueError(msg)


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def fnum(x):
    """float, or None for NaN / missing (JSON has no NaN)."""
    return None if x is None or not np.isfinite(float(x)) else float(x)


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False, default=str))


def fresh(path, is_dir=True):
    path = Path(path)
    check(not path.exists(), f'Output exists: {path}')
    if is_dir:
        path.mkdir(parents=True)
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
    return path


# ------------------------------------------------------------------ inputs
def read_input(folder):
    folder = Path(folder)
    lr = pd.read_csv(folder/'lr.csv', dtype=str)
    ent = pd.read_csv(folder/'entities.csv', dtype=str)
    check({'ligand', 'receptor', 'type'} <= set(lr) and not lr.duplicated(['ligand', 'receptor']).any(),
          f'{folder}/lr.csv: needs unique ligand,receptor,type rows')
    check({'entity', 'subunits'} <= set(ent) and ent.entity.is_unique, f'{folder}/entities.csv: entity,subunits')
    subunits = {r.entity: r.subunits.split(';') for r in ent.itertuples(index=False)}
    missing = sorted((set(lr.ligand) | set(lr.receptor)) - set(subunits))
    check(not missing, f'{folder}: entities without subunits {missing[:5]}')
    return lr, subunits, {'lr.csv': digest(folder/'lr.csv'), 'entities.csv': digest(folder/'entities.csv')}


def panel(inputs):
    """Union of the LR pairs of all inputs; entity definitions must agree."""
    lrs, subunits, hashes = [], {}, {}
    for folder in inputs:
        lr, sub, h = read_input(folder)
        for e, s in sub.items():
            check(subunits.get(e, s) == s, f'Entity {e} has different subunits in different inputs')
            subunits[e] = s
        lrs.append(lr[['ligand', 'receptor', 'type']])
        hashes[str(Path(folder).resolve())] = h
    lr = pd.concat(lrs).drop_duplicates()
    check(not lr.duplicated(['ligand', 'receptor']).any(), 'A pair has different types in different inputs')
    return lr.sort_values(['ligand', 'receptor']).reset_index(drop=True), subunits, hashes


def genes(args):
    lr, subunits, hashes = panel(args.input)
    out = fresh(args.output, is_dir=False)
    rows = {}
    for side in ('ligand', 'receptor'):
        for e in sorted(set(lr[side])):
            for g in subunits[e]:
                r = rows.setdefault(g, dict(symbol=g, as_ligand=False, as_receptor=False, entities=set()))
                r['as_' + side] = True
                r['entities'].add(e)
    table = pd.DataFrame([dict(r, entities=';'.join(sorted(r['entities']))) for r in rows.values()])
    table.sort_values('symbol').to_csv(out, index=False)
    write_json(str(out) + '.json', dict(created_utc=now(), inputs=hashes, n_pairs=len(lr),
               n_entities=len(set(lr.ligand) | set(lr.receptor)), n_genes=len(table),
               n_complex_entities=sum(len(subunits[e]) > 1 for e in set(lr.ligand) | set(lr.receptor))))
    print(f'{len(table)} subunit genes from {len(lr)} pairs -> {out}')


# ------------------------------------------------------------------ proteins
def read_fasta(path):
    seqs, name = {}, None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].strip(); seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: ''.join(v) for k, v in seqs.items()}


def load_proteins(folder):
    folder = Path(folder)
    meta = json.loads((folder/'proteins.json').read_text())
    check(meta.get('status') == 'complete' and digest(folder/'proteins.fasta') == meta['fasta_sha256'],
          'Protein map incomplete or changed')
    pmap = pd.read_csv(folder/'protein_map.csv', dtype=str)
    fasta = {k.split('|')[1]: v for k, v in read_fasta(folder/'proteins.fasta').items()}
    acc = dict(zip(pmap.lr_symbol, pmap.accession))
    check(set(acc) == set(fasta), 'protein_map.csv and proteins.fasta list different genes')
    return acc, fasta, meta


# ------------------------------------------------------------------ jobs
def subunit_pairs(ligand, receptor, subunits):
    """Ligand-subunit x receptor-subunit combinations, a protein with itself skipped. A pair left with no
    combination (homophilic, e.g. Cdh1-Cdh1) is scored as the protein with itself (reported separately)."""
    sp = [(a, b) for a in subunits[ligand] for b in subunits[receptor] if a != b]
    return sp or [(a, a) for a in subunits[ligand] if a in subunits[receptor]]


def pair_key(acc_a, acc_b):
    return '__'.join(sorted((acc_a, acc_b)))


def entity_table(lr, subunits, acc):
    rows = []
    for r in lr.itertuples(index=False):
        sp = subunit_pairs(r.ligand, r.receptor, subunits)
        unmapped = sorted({g for p in sp for g in p if g not in acc})
        keys = sorted({pair_key(acc[a], acc[b]) for a, b in sp if a in acc and b in acc})
        rows.append(dict(ligand=r.ligand, receptor=r.receptor, type=r.type,
                         homophilic=all(a == b for a, b in sp) and bool(sp),
                         n_subunit_pairs=len(sp), unmapped_subunits=';'.join(unmapped),
                         complete=bool(sp) and not unmapped and len(keys) > 0,
                         subunit_pair_keys=';'.join(keys)))
    return pd.DataFrame(rows)


def annotated_pairs(cellchat=None, exclude_csv=None):
    """(ligand, receptor) entity pairs annotated anywhere: CellChatDB (all classes) and/or a CSV."""
    pairs = set()
    if cellchat:
        from cellchat_inventory import load_cellchat
        it, _ = load_cellchat(cellchat)
        pairs |= set(zip(it.ligand, it.receptor))
    if exclude_csv:
        t = pd.read_csv(exclude_csv, dtype=str)
        pairs |= set(zip(t.ligand, t.receptor))
    return pairs


def cross_pairs(lr, subunits, annotated, n, seed):
    """n random unannotated combinations: panel ligand x panel receptor of the same type, not in the panel,
    not annotated, not the same entity, no shared subunit. Sampled without replacement."""
    panel_pairs = set(zip(lr.ligand, lr.receptor))
    cand = []
    for t, g in lr.groupby('type'):
        for lig in sorted(set(g.ligand)):
            for rec in sorted(set(g.receptor)):
                if (lig, rec) in panel_pairs or (lig, rec) in annotated or lig == rec:
                    continue
                if set(subunits[lig]) & set(subunits[rec]):
                    continue
                cand.append((lig, rec, t))
    check(len(cand) >= n, f'Only {len(cand)} unannotated combinations for {n} requested')
    rng = np.random.default_rng(seed)
    pick = sorted(rng.choice(len(cand), size=n, replace=False))
    return pd.DataFrame([cand[i] for i in pick], columns=['ligand', 'receptor', 'type']), len(cand)


def jobs(args):
    lr, subunits, hashes = panel(args.input)
    acc, fasta, pmeta = load_proteins(args.proteins)
    ent = entity_table(lr, subunits, acc).assign(set='cognate')
    if args.limit_pairs:
        rng = np.random.default_rng(args.seed)
        ok = ent[ent.complete]
        pick = set(rng.choice(ok.index, size=min(args.limit_pairs, len(ok)), replace=False))
        if args.include_longest:
            longest = max(ok.index, key=lambda i: max(
                sum(len(fasta[g]) for g in (a, b)) for a, b in subunit_pairs(ent.ligand[i], ent.receptor[i], subunits)))
            pick.add(longest)
        ent = ent.loc[sorted(pick)]
    n_cand = None
    ratio, cellchat = getattr(args, 'cross_ratio', 0.0), getattr(args, 'cellchat', None)
    exclude, cseed = getattr(args, 'exclude_pairs', None), getattr(args, 'cross_seed', 1)
    if ratio:
        check(cellchat or exclude, '--cross-ratio needs --cellchat (and/or --exclude-pairs)')
        n_cross = int(round(ratio * int(ent.complete.sum())))
        cr, n_cand = cross_pairs(lr, subunits, annotated_pairs(cellchat, exclude), n_cross, cseed)
        ent = pd.concat([ent, entity_table(cr, subunits, acc).assign(set='cross')], ignore_index=True)
    by_acc = {}
    for g, a in acc.items():
        check(by_acc.get(a, fasta[g]) == fasta[g], f'{a}: different sequences for one accession')
        by_acc[a] = fasta[g]
    keys = sorted({k for s in ent.subunit_pair_keys if s for k in s.split(';')})
    out = fresh(args.output)
    rows = []
    for k in keys:
        a, b = k.split('__')
        for order, (q, t) in zip(ORDERS, ((a, b), (b, a))):
            rows.append(dict(pair_id=f'{k}__{order}', key=k, order=order, query_accession=q, text_accession=t,
                             n_tokens=len(by_acc[q]) + len(by_acc[t]) + 3, query=by_acc[q], text=by_acc[t]))
    jb = pd.DataFrame(rows)
    jb.to_csv(out/'jobs.csv', index=False)
    jb[['pair_id', 'query', 'text']].to_csv(out/'plm_input.csv', index=False)
    ent.to_csv(out/'entity_pairs.csv', index=False)
    info = dict(created_utc=now(), inputs=hashes, proteins=str(Path(args.proteins).resolve()),
                proteins_fasta_sha256=pmeta['fasta_sha256'], unresolved_genes=pmeta.get('unresolved', []),
                limit_pairs=args.limit_pairs, seed=args.seed, include_longest=bool(args.include_longest),
                n_entity_pairs=int((ent.set == 'cognate').sum()),
                n_entity_pairs_complete=int((ent.complete & (ent.set == 'cognate')).sum()),
                cross_ratio=ratio, cross_seed=cseed, cellchat=cellchat, exclude_pairs=exclude,
                n_cross_pairs=int((ent.set == 'cross').sum()), n_cross_candidates=n_cand,
                n_cross_complete=int((ent.complete & (ent.set == 'cross')).sum()),
                n_subunit_pairs=len(keys), n_jobs=len(jb),
                max_tokens=int(jb.n_tokens.max()) if len(jb) else 0,
                n_jobs_over_official_max_tokens=int((jb.n_tokens > OFFICIAL_MAX_TOKENS).sum()),
                plm_input_sha256=digest(out/'plm_input.csv'))
    write_json(out/'jobs.json', info)
    print(json.dumps({k: info[k] for k in ['n_entity_pairs', 'n_entity_pairs_complete', 'n_cross_pairs',
                                           'n_cross_complete', 'n_subunit_pairs', 'n_jobs', 'max_tokens',
                                           'n_jobs_over_official_max_tokens']}, indent=2))


# ------------------------------------------------------------------ import
def auroc(pos, neg):
    pos, neg = np.asarray(pos, float), np.asarray(neg, float)
    if not len(pos) or not len(neg):
        return None
    r = pd.Series(np.concatenate([pos, neg])).rank().to_numpy()
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def separation(ppi, cross, n_boot=2000, seed=0):
    """AUROC panel pairs vs unannotated combinations (not negatives), with a bootstrap 95% interval."""
    rng = np.random.default_rng(seed)
    out = {}
    for name, (p, c) in {'all': (ppi, cross), **{t: (ppi[ppi.type == t], cross[cross.type == t])
                                                   for t in sorted(set(ppi.type))}}.items():
        a, b = p.score.dropna().to_numpy(), c.score.dropna().to_numpy()
        val = auroc(a, b)
        boots = [auroc(rng.choice(a, len(a)), rng.choice(b, len(b))) for _ in range(n_boot)] if val is not None else []
        out[name] = dict(n_panel=len(a), n_unannotated=len(b), auroc=val,
                         ci95=[fnum(np.quantile(boots, .025)), fnum(np.quantile(boots, .975))] if boots else None,
                         frac_gt_0_5_unannotated=fnum((b > 0.5).mean()) if len(b) else None)
    return out


def import_scores(args):
    jd = Path(args.jobs)
    info = json.loads((jd/'jobs.json').read_text())
    check(digest(jd/'plm_input.csv') == info['plm_input_sha256'], 'Jobs changed after export')
    meta = json.loads(Path(str(args.scores) + '.meta.json').read_text())
    check(meta.get('status') == 'complete' and meta.get('input_sha256') == info['plm_input_sha256'],
          'Scores were not produced from these jobs')
    if not args.exploratory:
        check(meta.get('checkpoint_sha256') == EXPECTED_CHECKPOINT_SHA256,
              'Checkpoint is not the registered humanV12 file (use --exploratory for others)')
    jb = pd.read_csv(jd/'jobs.csv', dtype=str)
    sc = pd.read_csv(args.scores, dtype={'pair_id': str})
    check(sc.pair_id.is_unique and set(sc.pair_id) == set(jb.pair_id), 'Every job must be scored exactly once')
    m = jb.merge(sc[['pair_id', 'score', 'n_tokens']], on='pair_id', suffixes=('', '_model'))
    check((m.n_tokens.astype(int) == m.n_tokens_model.astype(int)).all(), 'Token counts differ (truncated input?)')
    check(m.score.between(0, 1).all(), 'Scores outside [0,1]')
    w = m.pivot(index='key', columns='order', values='score')
    w['score'] = w[list(ORDERS)].mean(axis=1)
    w['abs_order_diff'] = (w.AB - w.BA).abs()
    ent = pd.read_csv(jd/'entity_pairs.csv', dtype={'unmapped_subunits': str, 'subunit_pair_keys': str})
    if 'set' not in ent:
        ent['set'] = 'cognate'
    rows = []
    for r in ent.itertuples(index=False):
        keys = r.subunit_pair_keys.split(';') if isinstance(r.subunit_pair_keys, str) and r.subunit_pair_keys else []
        s = w.loc[keys, 'score'] if keys else pd.Series(dtype=float)
        ok = bool(r.complete) and len(keys) > 0
        rows.append(dict(set=r.set, ligand=r.ligand, receptor=r.receptor, type=r.type, homophilic=bool(r.homophilic),
                         score=float(s.max()) if ok else np.nan,
                         best_subunit_pair=str(s.idxmax()) if ok else '',
                         n_scored_subunit_pairs=len(keys), n_subunit_pairs=r.n_subunit_pairs,
                         complete=ok, unmapped_subunits=r.unmapped_subunits if isinstance(r.unmapped_subunits, str) else ''))
    allp = pd.DataFrame(rows)
    ppi, cross = allp[allp.set == 'cognate'].drop(columns='set'), allp[allp.set == 'cross'].drop(columns='set')
    out = fresh(args.output)
    ppi.to_csv(out/'ppi_scores.csv', index=False)
    if len(cross):
        cross.to_csv(out/'cross_scores.csv', index=False)
    w.reset_index().to_csv(out/'subunit_pair_scores.csv', index=False)
    s = ppi.score.dropna()
    keep = ['checkpoint_name', 'checkpoint_revision', 'checkpoint_sha256', 'base_model', 'torch', 'transformers',
            'device', 'started_utc', 'finished_utc']
    report = dict(created_utc=now(), status='exploratory' if args.exploratory else 'registered-checkpoint',
                  jobs=info, scores_sha256=digest(args.scores), ppi_sha256=digest(out/'ppi_scores.csv'),
                  model={k: meta.get(k) for k in keep},
                  n_pairs=len(ppi), n_scored=int(ppi.complete.sum()),
                  n_missing=int((~ppi.complete).sum()), n_homophilic=int(ppi.homophilic.sum()),
                  homophilic_median=fnum(ppi.score[ppi.homophilic].median()) if ppi.homophilic.any() else None,
                  by_type={t: dict(n=int(len(d)), scored=int(d.complete.sum()),
                                   median=fnum(d.score.median()) if d.complete.any() else None,
                                   frac_gt_0_5=fnum((d.score.dropna() > 0.5).mean()) if d.complete.any() else None)
                           for t, d in ppi.groupby('type')},
                  score_distribution=dict(q05_25_50_75_95=[float(x) for x in s.quantile([.05, .25, .5, .75, .95])],
                                          frac_gt_0_99=float((s > 0.99).mean()), frac_lt_0_01=float((s < 0.01).mean()),
                                          frac_gt_0_5=float((s > 0.5).mean())) if len(s) else None,
                  order_asymmetry=dict(mean=float(w.abs_order_diff.mean()), max=float(w.abs_order_diff.max())),
                  panel_vs_unannotated=separation(ppi, cross) if len(cross) else None,
                  note='Model predictions, not measured affinities. Missing = some subunit unmapped; never imputed.')
    write_json(out/'score_report.json', report)
    print(json.dumps({k: report[k] for k in ['n_pairs', 'n_scored', 'n_missing', 'by_type', 'score_distribution',
                                            'panel_vs_unannotated']}, indent=2))


# ------------------------------------------------------------------ designs
def normalise(s, groups):
    """s / sum(s) within each group; NaN for a whole group if any member is NaN or the sum is 0."""
    s = pd.Series(s, dtype=float)
    tot = s.groupby(groups).transform(lambda x: x.sum() if x.notna().all() else np.nan)
    return (s / tot.where(tot > 0)).to_numpy()


def design_table(lr, ppi):
    d = lr[['ligand', 'receptor', 'type']].merge(ppi[['ligand', 'receptor', 'score']], on=['ligand', 'receptor'],
                                                 how='left', validate='one_to_one')
    d['k_ligand'] = d.groupby('ligand').receptor.transform('size')
    d['k_receptor'] = d.groupby('receptor').ligand.transform('size')
    d['s_raw'] = d.score
    d['s_norm_ligand'] = normalise(d.score, d.ligand)
    d['s_uniform_ligand'] = 1.0 / d.k_ligand
    d['s_norm_receptor'] = normalise(d.score, d.receptor)
    return d


def describe(d):
    ok = d.s_norm_ligand.notna()
    multi = d[ok & (d.k_ligand > 1)]
    spread = multi.groupby('ligand').s_raw.agg(lambda x: x.max() - x.min())
    y, x = -np.log(multi.s_norm_ligand), np.log(multi.k_ligand)
    r2 = float(np.corrcoef(x, y)[0, 1] ** 2) if len(multi) > 2 and x.std() > 0 and y.std() > 0 else None
    return dict(n_pairs=len(d), n_with_score=int(d.score.notna().sum()), n_with_s_norm_ligand=int(ok.sum()),
                n_ligands=int(d.ligand.nunique()), ligands_k1=int((d.groupby('ligand').size() == 1).sum()),
                pairs_in_k1_ligands=int((d.k_ligand == 1).sum()),
                ligands_with_missing_score=int(d.groupby('ligand').score.apply(lambda x: x.isna().any()).sum()),
                k_ligand_median=float(d.drop_duplicates('ligand').k_ligand.median()),
                k_ligand_max=int(d.k_ligand.max()),
                multi_receptor_ligands=int(multi.ligand.nunique()),
                within_ligand_raw_spread_lt_0_01=int((spread < 0.01).sum()),
                within_ligand_raw_spread_median=fnum(spread.median()) if len(spread) else None,
                r2_neglog_snorm_vs_log_k=r2,
                neglog_snorm_ligand_q05_50_95=[fnum(v) for v in (-np.log(d.s_norm_ligand[ok])).quantile([.05, .5, .95])] if ok.any() else None)


def designs(args):
    lr, _, h = read_input(args.input)
    ppi = pd.read_csv(args.ppi)
    d = design_table(lr, ppi)
    out = fresh(args.output, is_dir=False)
    d.to_csv(out, index=False)
    summary = dict(created_utc=now(), input=str(Path(args.input).resolve()), input_sha256=h,
                   ppi=str(Path(args.ppi).resolve()), ppi_sha256=digest(args.ppi), output_sha256=digest(out),
                   overall=describe(d), by_type={t: describe(g) for t, g in d.groupby('type')})
    write_json(str(out) + '.json', summary)
    print(json.dumps(summary['overall'], indent=2))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    g = sub.add_parser('genes'); g.add_argument('--input', nargs='+', action='extend', required=True); g.add_argument('--output', required=True)
    j = sub.add_parser('jobs'); j.add_argument('--input', nargs='+', action='extend', required=True)
    j.add_argument('--proteins', required=True); j.add_argument('--output', required=True)
    j.add_argument('--limit-pairs', type=int, default=0, help='smoke test: this many entity pairs (0 = all)')
    j.add_argument('--include-longest', action='store_true'); j.add_argument('--seed', type=int, default=0)
    j.add_argument('--cross-ratio', type=float, default=0.0, help='unannotated combinations per panel pair (plan: 3)')
    j.add_argument('--cross-seed', type=int, default=1)
    j.add_argument('--cellchat', help='CellChatDB folder (cellchat_inventory format) for the exclusion')
    j.add_argument('--exclude-pairs', help='extra CSV ligand,receptor to exclude (tests)')
    i = sub.add_parser('import'); i.add_argument('--jobs', required=True); i.add_argument('--scores', required=True)
    i.add_argument('--output', required=True); i.add_argument('--exploratory', action='store_true')
    d = sub.add_parser('designs'); d.add_argument('--input', required=True, help='one slice (lr.csv)')
    d.add_argument('--ppi', required=True); d.add_argument('--output', required=True)
    return p


if __name__ == '__main__':
    a = parser().parse_args()
    dict(genes=genes, jobs=jobs, designs=designs)[a.command](a) if a.command != 'import' else import_scores(a)
