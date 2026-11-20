"""Experiment 4 (reports/bindpred_ot_01/plan.md): BindPred affinities as S' = S / sum_R S in the OT cost.

  jobs       stage-2 inputs + protein map + frozen sequences -> oriented subunit pairs (ligand first)
  import     BindPred predictions -> ppi-style tables (pair score = best subunit pair) for stage2_ppi.py designs
  cases      case list per slice (core, permutations, sensitivity)
  check      reproducibility gate: no_ppi / V12 S' (beta 1) equal to the committed experiment-2 values
  summarise  decisions F1-F3 (plan section 4)

The OT runs use `stage2_s.py run --setting main` (eps 0.2, rho 1, 50 000 iterations, CPU) unchanged; this file
only builds its inputs and reads its outputs. BindPred scores are produced by `affinity.py embed/predict`.
S = 10 ** (pKd_hat - 10), i.e. 1/Kd in units of 1e10 per M; only the per-ligand normalised S' enters the OT, so the
unit cancels.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import pandas as pd

from research import check

N_PERM = 100                          # plan section 3
ALPHA = 0.05
SCALE_PKD = 10.0                      # S = 10 ** (pKd_hat - SCALE_PKD); cancels in S'
CORE = ['no_ppi', 'bp_snorm', 'uniform']
GATE_TOL = 1e-10                      # plan section 6: reproduction of experiment-2 values


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False,
                                     default=lambda x: x.item() if hasattr(x, 'item') else str(x)))


def digest(path):
    import hashlib
    h = hashlib.sha256()
    with Path(path).open('rb') as s:
        for block in iter(lambda: s.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# ------------------------------------------------------------------ jobs
def jobs(args):
    """Oriented (ligand subunit, receptor subunit) pairs of the union panel, same combinations as stage2_ppi.py."""
    import stage2_ppi as sp
    out = Path(args.output_dir)
    check(not out.exists(), f'{out} exists')
    lr, subunits, hashes = sp.panel(args.input)
    pm = pd.read_csv(args.protein_map, dtype=str)
    acc = dict(zip(pm.lr_symbol, pm.accession))
    seq = pd.read_csv(args.sequences, sep='\t', dtype={'length': int})
    import hashlib
    check((seq.sequence.map(lambda s: hashlib.sha256(s.encode()).hexdigest()) == seq.sequence_sha256).all(),
          'Frozen sequences do not match their sha256')
    known = dict(zip(pm.drop_duplicates('accession').accession, pm.drop_duplicates('accession').sequence_sha256))
    seqs = dict(zip(seq.accession, seq.sequence))
    check(all(hashlib.sha256(seqs[a].encode()).hexdigest() == h for a, h in known.items() if a in seqs)
          and set(known) <= set(seqs), 'Sequences differ from the stage-2 protein map')
    ent, oriented = [], {}
    for r in lr.itertuples(index=False):
        combos = sp.subunit_pairs(r.ligand, r.receptor, subunits)
        unmapped = sorted({g for c in combos for g in c if g not in acc})
        keys = sorted({f'{acc[a]}__{acc[b]}' for a, b in combos if a in acc and b in acc})   # ligand first
        for k in keys:
            oriented[k] = k.split('__')
        ent.append(dict(ligand=r.ligand, receptor=r.receptor, type=r.type,
                        homophilic=all(a == b for a, b in combos) and bool(combos),
                        complete=bool(combos) and not unmapped and len(keys) > 0,
                        unmapped_subunits=';'.join(unmapped), oriented_keys=';'.join(keys)))
    ent = pd.DataFrame(ent)
    pairs = pd.DataFrame([dict(pair_id=k, lig_acc=a, rec_acc=b, lig_seq=seqs[a], rec_seq=seqs[b])
                          for k, (a, b) in sorted(oriented.items())])
    out.mkdir(parents=True)
    ent.to_csv(out/'entity_pairs.csv', index=False)
    pairs.to_csv(out/'pairs.csv', index=False)
    write_json(out/'jobs.json', dict(created_utc=now(), inputs=hashes, protein_map_sha256=digest(args.protein_map),
               sequences_sha256=digest(args.sequences), n_entity_pairs=len(ent), n_complete=int(ent.complete.sum()),
               n_oriented_pairs=len(pairs), n_sequences=int(len(set(pairs.lig_seq) | set(pairs.rec_seq))),
               n_longer_than_1022=int(sum(len(s) > 1022 for s in set(pairs.lig_seq) | set(pairs.rec_seq)))))
    print(f'{len(ent)} entity pairs ({int(ent.complete.sum())} complete), {len(pairs)} oriented subunit pairs -> {out}')


# ------------------------------------------------------------------ import
def import_scores(args):
    """Pair score = best (highest pKd_hat) subunit pair, as for V12 (max over subunit pairs)."""
    out = Path(args.output_dir)
    check(not out.exists(), f'{out} exists')
    ent = pd.read_csv(Path(args.jobs_dir)/'entity_pairs.csv', dtype={'oriented_keys': str, 'unmapped_subunits': str})
    bp = pd.read_csv(args.bindpred).set_index('pair_id').pkd_hat
    rows = []
    for r in ent.itertuples(index=False):
        keys = r.oriented_keys.split(';') if isinstance(r.oriented_keys, str) and r.oriented_keys else []
        check(all(k in bp.index for k in keys), f'{r.ligand}-{r.receptor}: subunit pair without prediction')
        ok = bool(r.complete)
        best = bp[keys].idxmax() if ok else ''
        pk = float(bp[keys].max()) if ok else np.nan
        sorted_key = '__'.join(sorted(best.split('__'))) if ok else ''          # stage2_s.pair_table format
        rows.append(dict(ligand=r.ligand, receptor=r.receptor, type=r.type, homophilic=bool(r.homophilic),
                         pkd_hat=pk, score=10 ** (pk - SCALE_PKD) if ok else np.nan,
                         best_subunit_pair=sorted_key, best_oriented=best, complete=ok))
    ppi = pd.DataFrame(rows)
    # sensitivity (plan section 5): receptor-only score = mean pKd_hat of the receptor over its scored panel pairs
    rec_mean = ppi.groupby('receptor').pkd_hat.mean()
    rec = ppi.assign(pkd_hat=ppi.receptor.map(rec_mean).where(ppi.complete))
    rec['score'] = 10 ** (rec.pkd_hat - SCALE_PKD)
    out.mkdir(parents=True)
    ppi.to_csv(out/'ppi_bindpred.csv', index=False)
    rec.to_csv(out/'ppi_bindpred_receptor_only.csv', index=False)
    multi = ppi[ppi.complete].groupby('ligand').filter(lambda g: len(g) > 1)
    spread = multi.groupby('ligand').pkd_hat.agg(lambda x: x.max() - x.min())
    write_json(out/'import.json', dict(created_utc=now(), bindpred_sha256=digest(args.bindpred), n_pairs=len(ppi),
               n_scored=int(ppi.complete.sum()), pkd_hat=ppi.pkd_hat.describe().to_dict(),
               multi_receptor_ligands=int(spread.size), within_ligand_pkd_spread_median=float(spread.median()),
               within_ligand_pkd_spread_lt_0_1=int((spread < 0.1).sum())))
    print(f'{int(ppi.complete.sum())}/{len(ppi)} pairs scored -> {out}')


# ------------------------------------------------------------------ cases
def case(name, s, beta, perm=0, group='core', keep=True):
    return dict(name=name, s=s, beta=beta, perm=perm, override=None, panel='all', nitermax=50000, keep=keep, group=group)


def cases_cmd(args):
    """Plan section 3: one case list per (slice, designs). Names are unique across the three lists."""
    n_perm = 2 if args.smoke else N_PERM
    lists = {'bp': [case('no_ppi', 'none', 0.0), case('bp_snorm', 'snorm', 1.0), case('uniform', 'uniform', 1.0)]
                   + [case(f'bp_snorm_perm_{k:03d}', 'snorm', 1.0, perm=k, group='perm', keep=False)
                      for k in range(1, n_perm + 1)],
             'bprec': [case('bp_receptor_snorm', 'snorm', 1.0, group='sensitivity')],
             'v12': [case('v12_snorm', 'snorm', 1.0, group='sensitivity')]}
    out = Path(args.runs)
    for sl in ('post', 'ant'):
        for tag, cs in lists.items():
            f = out/f'{sl}_{tag}'/'cases.json'
            check(not f.exists(), f'{f} exists')
            f.parent.mkdir(parents=True, exist_ok=True)
            write_json(f, dict(created_utc=now(), setting='main', slice=sl, smoke=args.smoke, cases=cs))
    print(f'case lists -> {out} ({len(lists["bp"])} + 1 + 1 per slice)')


# ------------------------------------------------------------------ reproducibility gate
def case_values(folder, name, top_tf=10):
    """Receptor partial Spearman of one finished case (the rest of the folder may not have run yet)."""
    rec = json.loads((Path(folder)/name/'case.json').read_text())
    check(not rec['error'], f'{folder}/{name}: solver error {rec["error"]}')
    m = pd.read_csv(Path(folder)/name/'metrics.csv')
    return m[m.top_tf.astype(str) == str(top_tf)].set_index('receptor').partial_spearman


def check_cmd(args):
    """no_ppi and V12 S' (beta 1) must reproduce experiment 2 (columns no_ppi, snorm_b1) to GATE_TOL."""
    import stage2_s as ss
    res = []
    for sl in ('post', 'ant'):
        ref = pd.read_csv(Path(args.exp2_summary)/f'{sl}_main_receptor_partial_spearman.csv', index_col=0)
        for folder, ours, theirs in [(f'{sl}_bp', 'no_ppi', 'no_ppi'), (f'{sl}_v12', 'v12_snorm', 'snorm_b1')]:
            a, b = case_values(Path(args.runs)/folder, ours), ref[theirs]
            check(set(a.index) == set(b.index), f'{sl}/{ours}: receptor sets differ from experiment 2')
            b = b.loc[a.index]
            nan_mismatch = int((a.isna() != b.isna()).sum())
            diff = float((a - b).abs().max()) if (a.notna() & b.notna()).any() else 0.0
            res.append(dict(slice=sl, case=ours, reference=theirs, receptors=int(len(a)),
                            max_abs_diff=diff, nan_mismatch=nan_mismatch))
    passed = all(r['nan_mismatch'] == 0 and r['max_abs_diff'] <= GATE_TOL for r in res)
    write_json(args.output, dict(created_utc=now(), tolerance=GATE_TOL, cases=res, passed=passed))
    print(json.dumps(res), 'passed' if passed else 'FAILED')
    check(passed, 'Experiment-2 values were not reproduced')


# ------------------------------------------------------------------ summarise
def slice_result(runs, sl, top_tf):
    import stage2_s as ss
    wide, info = ss.load_runs(Path(runs)/f'{sl}_bp', top_tf)
    for tag in ('bprec', 'v12'):
        w, i = ss.load_runs(Path(runs)/f'{sl}_{tag}', top_tf)
        wide = wide.join(w, how='outer'); info.update(i)
    core = [n for n, c in info.items() if c['group'] in ('core', 'perm')]
    failed = [n for n in core if info[n]['error']]
    common = wide[[n for n in core if n in wide]].dropna()
    M = common.mean()
    perms = sorted(n for n in M.index if n.startswith('bp_snorm_perm_'))
    mp = M[perms].to_numpy()
    out = dict(receptors_common=int(len(common)), failed_cases=failed,
               undetermined=any(info[n]['error'] for n in CORE if n in info),
               M={k: float(M[k]) for k in CORE if k in M})
    if 'bp_snorm' in M and len(mp):
        p = float((1 + (mp >= M['bp_snorm'] - 1e-12).sum()) / (1 + len(mp)))
        out['F1'] = dict(M=float(M['bp_snorm']), n_perm=len(mp), perm_mean=float(mp.mean()),
                         perm_sd=float(mp.std(ddof=1)) if len(mp) > 1 else None, p=p,
                         label='information' if p <= ALPHA else 'no_information')
    out['F2_vs_uniform'] = ss.paired(common, 'bp_snorm', 'uniform')
    out['F3_vs_no_ppi'] = ss.paired(common, 'bp_snorm', 'no_ppi')
    out['uniform_vs_no_ppi'] = ss.paired(common, 'uniform', 'no_ppi')
    sens = wide[[c for c in ['bp_snorm', 'uniform', 'no_ppi', 'bp_receptor_snorm', 'v12_snorm'] if c in wide]]
    out['sensitivity'] = {
        'bp_receptor_snorm_vs_uniform': ss.paired(sens, 'bp_receptor_snorm', 'uniform'),
        'bp_snorm_vs_bp_receptor_snorm': ss.paired(sens, 'bp_snorm', 'bp_receptor_snorm'),
        'v12_snorm_vs_uniform': ss.paired(sens, 'v12_snorm', 'uniform'),
        'bp_snorm_vs_v12_snorm': ss.paired(sens, 'bp_snorm', 'v12_snorm')}
    return wide, out


def conclusion(f1, f2):
    if 'undetermined' in (f1, f2):
        return '无法判定(核心条件求解失败)'
    if f1 != 'information':
        return '看不出 BindPred 的配对信息'
    if f2 == 'improved':
        return 'BindPred 的 S′ 改进了推断,且超过平均分配'
    return '带有信息,但不超过平均分配'


def summarise(args):
    import stage2_s as ss
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    out.mkdir(parents=True)
    details = {}
    for sl in ('post', 'ant'):
        wide, res = slice_result(args.runs, sl, ss.TOP_TF[0])
        wide.to_csv(out/f'{sl}_receptor_partial_spearman.csv')
        for k in ss.TOP_TF[1:]:
            w, r = slice_result(args.runs, sl, k)
            res[f'top_tf_{k}'] = {key: r[key] for key in ('F2_vs_uniform', 'F3_vs_no_ppi')}
        details[sl] = res
    lab = lambda sl, key: ('undetermined' if details[sl]['undetermined'] else
                           (details[sl].get(key) or {}).get('label'))
    f1p, f1a, f2p, f2a = lab('post', 'F1'), lab('ant', 'F1'), lab('post', 'F2_vs_uniform'), lab('ant', 'F2_vs_uniform')
    verdict = dict(F1_post=f1p, F1_ant=f1a, F1_confirmed=bool(f1p == f1a),
                   F2_post=f2p, F2_ant=f2a, F2_confirmed=bool(f2p == f2a),
                   F3_post=lab('post', 'F3_vs_no_ppi'), F3_ant=lab('ant', 'F3_vs_no_ppi'),
                   conclusion_post=conclusion(f1p, f2p), conclusion_ant=conclusion(f1a, f2a),
                   failed_cases=details['post']['failed_cases'] + details['ant']['failed_cases'])
    write_json(out/'decisions.json', dict(created_utc=now(), n_perm=N_PERM, alpha=ALPHA, verdict=verdict,
                                          details=details))
    print(json.dumps(verdict, indent=2, ensure_ascii=False))


def parser():
    pa = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = pa.add_subparsers(dest='command', required=True)
    j = sub.add_parser('jobs')
    j.add_argument('--input', nargs='+', action='extend', required=True, help='stage-2 input folders')
    j.add_argument('--protein-map', required=True); j.add_argument('--sequences', required=True)
    j.add_argument('--output-dir', required=True)
    i = sub.add_parser('import')
    i.add_argument('--jobs-dir', required=True); i.add_argument('--bindpred', required=True)
    i.add_argument('--output-dir', required=True)
    c = sub.add_parser('cases'); c.add_argument('--runs', required=True); c.add_argument('--smoke', action='store_true')
    k = sub.add_parser('check'); k.add_argument('--runs', required=True); k.add_argument('--exp2-summary', required=True)
    k.add_argument('--output', required=True)
    s = sub.add_parser('summarise'); s.add_argument('--runs', required=True); s.add_argument('--output', required=True)
    return pa


if __name__ == '__main__':
    a = parser().parse_args()
    dict(jobs=jobs, cases=cases_cmd, check=check_cmd, summarise=summarise)[a.command](a) \
        if a.command != 'import' else import_scores(a)
