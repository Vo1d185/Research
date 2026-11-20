"""Experiment 2 (reports/stage2_ppi_02/plan.md): real slices, no S vs raw V12 S vs the advisor's S' = S / sum_R S.

  cases        list the cases of one slice and setting (conditions, permutations, sensitivity) -> cases.json
  run          solve cases (by name) and evaluate them at once: receptor-level received signal -> partial
               Spearman against the scSeqComm target score (top 10 / 5 / 20 TFs, mean), with the functions of
               evaluate_targets.py. Non-permutation cases also keep the signal; --save-case writes a stage-2
               result folder (full transport) for the equivalence gate
  equivalence  compare this script's receptor metrics with evaluate_targets.py run on the saved folders
  summarise    D1 / D2 / D3 per slice and setting (posterior decides, anterior confirms), sensitivity, stability

Cost: C = -log F(d) - beta * log S_eff, F(d) = r0 / (r0 + d), r0 = 50 um, Stage 2 condition A otherwise.
Environments: commot-research (setting main, CPU); setting second on GPU: the environment with torch
(inject-gpu) and --device cuda (gpu_unot.py, amendment_01 of stage2_inject_01). All paths are arguments.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from research import check

R0_UM = 50.0
CUTOFF_UM = 363.0757               # plan section 2 (Stage 2 condition A); fixed here, not read from config.json (amendment 01)
STOPTHR = 1e-11
TOP_TF = (10, 5, 20)                 # first = primary rule (stage2_plan.md 9.3)
MAX_TOKENS = 1603
N_BOOT, BOOT_SEED, ALPHA = 1000, 1, 0.05
STABILITY_NITER = 200000
SETTINGS = {   # plan section 2; beta chosen on the injection experiment before any real evaluation
    'main': dict(eps=0.2, rho=1.0, nitermax=50000, beta_raw=1.0, beta_snorm=0.2, n_perm=100, device='cpu'),
    'second': dict(eps=0.1, rho=10.0, nitermax=50000, beta_raw=1.0, beta_snorm=1.0, n_perm=0, device='cuda'),
}
CORE = ['no_ppi', 'raw', 'snorm', 'uniform']


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


def load_input(source):
    """stage2_inject.load, with the diffusion edges built at the plan's cutoff (CUTOFF_UM) whatever config.json stores.
    stage2_input_post stores 500 (Stage 2 passed 363.0757 on the command line); amendment 01."""
    import stage2
    import stage2_inject as si
    prep, cfg, x, spatial, lr, ligands, receptors, a, b, edges = si.load(source)
    stored = float(cfg['diffusion_cutoff_um'])
    if stored != CUTOFF_UM:
        edges = dict(zip(['diffusion', 'contact'], stage2.spatial_edges(spatial, CUTOFF_UM)))
    cfg = dict(cfg, diffusion_cutoff_um=CUTOFF_UM, config_cutoff_um=stored)
    return prep, cfg, x, spatial, lr, ligands, receptors, a, b, edges


# ------------------------------------------------------------------ cases and scores
def case_list(setting, slice_, smoke=False):
    """Plan sections 3 and 5. Each case: name, s (none|raw|snorm|uniform), beta, perm (seed or 0),
    override (None|long|homo), panel (all|scored), nitermax, keep (save the signal)."""
    st = SETTINGS[setting]
    c = lambda name, s, beta, perm=0, override=None, panel='all', nitermax=st['nitermax'], keep=True, group='core': dict(
        name=name, s=s, beta=beta, perm=perm, override=override, panel=panel, nitermax=nitermax, keep=keep, group=group)
    out = [c('no_ppi', 'none', 0.0), c('raw', 'raw', st['beta_raw']), c('snorm', 'snorm', st['beta_snorm']),
           c('uniform', 'uniform', st['beta_snorm'])]
    n_perm = (2 if smoke else st['n_perm'])
    for k in range(1, n_perm + 1):
        out += [c(f'raw_perm_{k:03d}', 'raw', st['beta_raw'], perm=k, keep=False, group='perm'),
                c(f'snorm_perm_{k:03d}', 'snorm', st['beta_snorm'], perm=k, keep=False, group='perm')]
    if setting == 'main':
        out.append(c('snorm_b1', 'snorm', 1.0, group='sensitivity'))
        sens = [('long_median', 'long', 'all'), ('homo_median', 'homo', 'all'), ('scored_only', None, 'scored')]
        for tag, override, panel in (sens[:1] if smoke else sens):
            out += [c(f'raw_{tag}', 'raw', st['beta_raw'], override=override, panel=panel, group='sensitivity'),
                    c(f'snorm_{tag}', 'snorm', st['beta_snorm'], override=override, panel=panel, group='sensitivity')]
            if panel == 'scored':
                out.append(c('no_ppi_scored_only', 'none', 0.0, panel='scored', group='sensitivity'))
    if setting == 'second' and slice_ == 'post':
        out += [c('raw_n200k', 'raw', st['beta_raw'], nitermax=STABILITY_NITER, group='stability'),
                c('snorm_n200k', 'snorm', st['beta_snorm'], nitermax=STABILITY_NITER, group='stability')]
    if smoke and setting == 'second':
        out = [dict(x, nitermax=1000) for x in out if x['name'] == 'no_ppi']
    return out


def pair_table(lr, designs, ppi, protein_map):
    """lr order; raw V12 score (NaN if unscored), homophilic flag, tokens of the best subunit pair."""
    d = lr[['ligand', 'receptor', 'type']].merge(designs[['ligand', 'receptor', 's_raw']], on=['ligand', 'receptor'],
                                                 how='left', validate='one_to_one')
    check(len(d) == len(lr) and len(designs) == len(lr), 'Designs and input pair lists differ')
    d = d.merge(ppi[['ligand', 'receptor', 'homophilic', 'best_subunit_pair']], on=['ligand', 'receptor'], how='left')
    length = protein_map.drop_duplicates('accession').set_index('accession').length
    def tokens(key):
        if not isinstance(key, str) or '__' not in key:
            return np.nan
        a, b = key.split('__')
        return float(length.get(a, np.nan) + length.get(b, np.nan) + 3)
    d['tokens'] = d.best_subunit_pair.map(tokens)
    d['homophilic'] = d.homophilic.fillna(False).astype(bool)
    check(np.isfinite(d.tokens[d.s_raw.notna()]).all(), 'A scored pair has no sequence length')
    check((d.s_raw.dropna() > 0).all(), 'Scores must be positive')
    return d


def effective(pairs, spec):
    """S_eff per pair of the case panel (plan sections 3 and 5) and the panel mask over the input pairs."""
    raw = pairs.s_raw.to_numpy(float).copy()
    scored = np.isfinite(raw)
    median = {t: float(np.median(raw[scored & (pairs.type == t).to_numpy()])) for t in sorted(pairs.type.unique())}
    if spec['perm']:
        rng = np.random.default_rng(spec['perm'])
        for t in sorted(median):
            idx = np.flatnonzero(scored & (pairs.type == t).to_numpy())
            raw[idx] = raw[idx][rng.permutation(len(idx))]
    fill = pairs.type.map(median).to_numpy(float)
    s = np.where(scored, raw, fill)                               # unscored pairs: type median, never permuted
    if spec['override'] == 'long':
        s = np.where(scored & (pairs.tokens.to_numpy(float) > MAX_TOKENS), fill, s)
    elif spec['override'] == 'homo':
        s = np.where(scored & pairs.homophilic.to_numpy(bool), fill, s)
    panel = scored if spec['panel'] == 'scored' else np.ones(len(pairs), bool)
    sub = pairs[panel].reset_index(drop=True)
    s = s[panel]
    if spec['s'] == 'snorm':
        s = s / pd.Series(s).groupby(sub.ligand.to_numpy()).transform('sum').to_numpy()
    elif spec['s'] == 'uniform':
        s = 1.0 / sub.groupby('ligand').receptor.transform('size').to_numpy(float)
    elif spec['s'] == 'none':
        s = np.ones(len(sub))
    return s, panel


# ------------------------------------------------------------------ evaluation (same functions as evaluate_targets)
class Evaluation:
    def __init__(self, inp, genes_dir, db, top_tf=TOP_TF):
        import evaluate_targets as et
        self.et = et
        inp, genes_dir = Path(inp), Path(genes_dir)
        export = json.loads((genes_dir/'gene_export.json').read_text())
        check(export.get('status') == 'complete' and digest(genes_dir/'gene_expression.npz') ==
              export['gene_expression_sha256'], 'Gene export incomplete or changed')
        with np.load(genes_dir/'gene_expression.npz', allow_pickle=False) as z:
            X, spots, genes = z['X'], z['spots'].astype(str), z['genes'].astype(str)
        spatial = pd.read_csv(inp/'spatial.csv', index_col=0)
        check(list(spots) == list(spatial.index.astype(str)), 'Gene export spots differ from --input')
        table = pd.read_csv(genes_dir/'genes.csv', dtype={'data_symbol': str, 'lr_symbol': str})
        symbol = dict(zip(table.lr_symbol.dropna(), table.data_symbol[table.lr_symbol.notna()]))
        detected = dict(zip(table.data_symbol, table.detected_fraction))
        expression = pd.read_csv(inp/'expression.csv', index_col=0)
        expression.index = expression.index.astype(str)
        expression = expression.loc[spots]
        receptor_tf, tf_tg, self.db_meta = et.load_db(db)
        panel = pd.read_csv(inp/'lr.csv', dtype=str)
        subunits = {}
        if (inp/'entities.csv').is_file():
            ent = pd.read_csv(inp/'entities.csv', dtype=str)
            subunits = {r.entity: r.subunits.split(';') for r in ent.itertuples(index=False)}
        members = lambda e: subunits.get(e, [e])
        panel_genes = {symbol[g] for e in set(panel.ligand) | set(panel.receptor) for g in members(e) if g in symbol}
        receptors = sorted(set(panel.receptor))
        check(all(m in symbol for r in receptors for m in members(r)), 'Receptors/subunits without data symbol')
        sets = pd.concat([et.target_sets(receptors, symbol, receptor_tf, tf_tg, detected, panel_genes, k, 5, 0.05,
                                         subunits, 'any') for k in top_tf]).drop_duplicates(['receptor', 'top_tf'])
        self.sets = sets
        self.targets = {(r.receptor, str(r.top_tf)): et.target_scores(X, genes, r.targets.split(';'), 'mean')
                        for r in sets[sets.status == 'evaluated'].itertuples(index=False)}
        self.receptor_expr = {r: expression[r].to_numpy(float) for r in receptors}
        self.spots = list(spots)

    def metrics(self, per_receptor):
        rows = []
        for (receptor, k), target in self.targets.items():
            signal = per_receptor.get(receptor)
            if signal is None:
                continue
            rows.append(dict(receptor=receptor, top_tf=k, total_received=float(signal.sum()),
                             spearman=self.et.spearman(signal, target),
                             partial_spearman=self.et.partial_spearman(signal, target, [self.receptor_expr[receptor]])))
        return pd.DataFrame(rows)


# ------------------------------------------------------------------ run
def solve(row, col, c, a, b, st, nitermax, device):
    if device == 'cpu':
        from research import solve_sparse
        return solve_sparse(row, col, c, a, b, st['eps'], st['rho'], nitermax, STOPTHR)
    from gpu_unot import solve_sparse_gpu
    return solve_sparse_gpu(row, col, c, a, b, st['eps'], st['rho'], nitermax, STOPTHR, device=device)


def save_case(folder, prep, inp, case, lr_case, row, col, c, T, d, p, a, b, x, ligands, receptors, info):
    """A stage-2 result folder (stage2.py run format) that evaluate_targets.py can read."""
    import plans
    import stage2
    folder = Path(folder)
    check(not folder.exists(), f'{folder} exists')
    folder.mkdir(parents=True)
    lr_case.to_csv(folder/'used_lr.csv', index=False)
    plans.save_sparse(folder/'case_000.npz', row, col, c, T, (len(ligands)*len(x), len(receptors)*len(x)), a, b,
                      x.index, ligands, receptors, D=d, pair=p, pair_ligand=lr_case.ligand.to_numpy(str),
                      pair_receptor=lr_case.receptor.to_numpy(str), pair_type=lr_case.type.to_numpy(str))
    pd.DataFrame([dict(case='case_000', family='r0_raw', kernel='r0', r0=R0_UM, r0_um=R0_UM,
                       score_mode='no_ppi' if case['s'] == 'none' else 'ppi', niter=info.get('niter'))]
                 ).to_csv(folder/'summary.csv', index=False)
    write_json(folder/'manifest.json', dict(status='complete', method=f'stage2_s.py ({case["name"]})', h=prep['neighbour_spacing_um'],
               shuffle=prep.get('shuffle'), input=str(Path(inp).resolve()), cutoff_um=CUTOFF_UM,
               input_sha256={f: stage2.digest(Path(inp)/f) for f in ['expression.csv', 'spatial.csv', 'lr.csv']}))


def run(args):
    import plans
    import stage2
    import stage2_inject as si
    st = SETTINGS[args.setting]
    device = args.device or st['device']
    cases = {c['name']: c for c in json.loads(Path(args.cases).read_text())['cases']}
    todo = [cases[n] for n in args.case]
    prep, cfg, x, spatial, lr, ligands, receptors, a, b, edges = load_input(args.input)
    meta = json.loads(Path(str(args.designs) + '.json').read_text())
    check(meta['input_sha256']['lr.csv'] == digest(Path(args.input)/'lr.csv'), 'Designs were made for another input')
    pairs = pair_table(lr, pd.read_csv(args.designs), pd.read_csv(args.ppi_scores), pd.read_csv(args.protein_map))
    ev = Evaluation(args.input, args.genes, args.db)
    n = len(x)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    for case in todo:
        target = out/case['name']
        if (target/'case.json').exists():
            continue                                               # finished earlier (resumable)
        target.mkdir(exist_ok=True)
        s, panel = effective(pairs, case)
        lr_case = lr[panel].reset_index(drop=True)
        row, col, c0, d, p = stage2.build_sparse_cost(lr_case, ligands, receptors, a, b, n, edges, R0_UM)
        c = c0 + case['beta'] * (-np.log(s))[p]
        start = time.time()
        error = ''
        try:
            T, info = solve(row, col, c, a, b, st, case['nitermax'], device)
        except ValueError as e:                                    # solver checks failed: not converged
            T, info, error = None, {}, str(e)
        seconds = time.time() - start
        rec = dict(case=case, setting=args.setting, slice=args.slice, device=device, input=str(Path(args.input).resolve()),
                   designs_sha256=digest(args.designs), cutoff_um=cfg['diffusion_cutoff_um'],
                   config_cutoff_um=cfg['config_cutoff_um'], entries=int(len(c)), seconds=round(seconds, 1), error=error,
                   niter=info.get('niter'), cap_reached=info.get('cap_reached'),
                   checks={k: v for k, v in info.items() if k.startswith(('mass', 'kkt', 'max', 'row', 'shape'))},
                   total_transport=float(T.sum()) if T is not None else None,
                   s_eff=dict(min=float(s.min()), max=float(s.max()), n_pairs=int(len(s))))
        if T is not None:
            plan = dict(sparse=True, ligands=np.asarray(ligands, str), receptors=np.asarray(receptors, str),
                        row=row, col=col, T=T)
            _, per_receptor = plans.received_signal(plan, lr_case, n)
            m = ev.metrics(per_receptor)
            m.to_csv(target/'metrics.csv', index=False)
            if case['keep']:
                names = sorted(per_receptor)
                np.savez_compressed(target/'signal.npz', receptors=np.asarray(names, str),
                                    signal=np.vstack([per_receptor[r] for r in names]))
            if args.save_case:
                save_case(Path(args.save_case)/case['name'], prep, args.input, case, lr_case, row, col, c, T, d, p,
                          a, b, x, ligands, receptors, info)
        write_json(target/'case.json', rec)
        print(f'{now()} {args.slice}/{args.setting}/{case["name"]}: {len(c):,} entries, {seconds:.0f}s, '
              f'{info.get("niter")} it' + (f' ERROR {error}' if error else ''), flush=True)


def cases_cmd(args):
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    out.parent.mkdir(parents=True, exist_ok=True)
    cs = case_list(args.setting, args.slice, args.smoke)
    write_json(out, dict(created_utc=now(), setting=args.setting, slice=args.slice, smoke=args.smoke,
                         parameters=SETTINGS[args.setting], r0_um=R0_UM, stopthr=STOPTHR, cases=cs))
    print(f'{len(cs)} cases -> {out}')


# ------------------------------------------------------------------ equivalence gate
def equivalence(args):
    ev = pd.read_csv(Path(args.eval_dir)/'case_metrics.csv')
    ev = ev[(ev.level == 'receptor') & (ev.top_tf.astype(str) == str(TOP_TF[0])) & (ev.target_score == 'mean')
            & (ev.rule == 'ours')]
    res = []
    for name in args.case:
        ours = pd.read_csv(Path(args.runs)/name/'metrics.csv')
        ours = ours[ours.top_tf.astype(str) == str(TOP_TF[0])].set_index('receptor').partial_spearman
        theirs = ev[ev.results.map(lambda r: Path(r).name) == name].set_index('receptor').partial_spearman
        check(len(theirs) and set(ours.index) == set(theirs.index), f'{name}: receptor sets differ')
        diff = (ours - theirs.loc[ours.index]).abs()
        both_nan = ours.isna() & theirs.loc[ours.index].isna()
        res.append(dict(case=name, receptors=int(len(ours)), max_abs_diff=float(diff[~both_nan].max()),
                        nan_mismatch=int((ours.isna() != theirs.loc[ours.index].isna()).sum())))
    passed = all(r['max_abs_diff'] <= 1e-12 and r['nan_mismatch'] == 0 for r in res)
    write_json(args.output, dict(created_utc=now(), cases=res, tolerance=1e-12, passed=passed))
    print(json.dumps(res), 'passed' if passed else 'FAILED')
    check(passed, 'This script and evaluate_targets.py disagree')


# ------------------------------------------------------------------ summarise
def load_runs(folder, top_tf=TOP_TF[0]):
    folder = Path(folder)
    cases = json.loads((folder/'cases.json').read_text())['cases']
    wide, info = {}, {}
    for c in cases:
        rec = json.loads((folder/c['name']/'case.json').read_text())
        info[c['name']] = dict(c, error=rec['error'], niter=rec['niter'], cap_reached=rec['cap_reached'],
                               total_transport=rec['total_transport'], seconds=rec['seconds'])
        if not rec['error']:
            m = pd.read_csv(folder/c['name']/'metrics.csv')
            wide[c['name']] = m[m.top_tf.astype(str) == str(top_tf)].set_index('receptor').partial_spearman
    return pd.DataFrame(wide), info


def paired(wide, a, b, rng_seed=BOOT_SEED, n_boot=N_BOOT):
    """Mean over receptors defined in both of (a - b), receptor bootstrap; decision label."""
    if a not in wide or b not in wide:
        return dict(delta=None, ci95=[None, None], receptors=0, label='not_run')
    d = (wide[a] - wide[b]).dropna().to_numpy()
    if not len(d):
        return dict(delta=None, ci95=[None, None], receptors=0, label='no_receptor')
    rng = np.random.default_rng(rng_seed)
    boots = d[rng.integers(len(d), size=(n_boot, len(d)))].mean(1)
    lo, hi = (float(x) for x in np.percentile(boots, [2.5, 97.5]))
    label = 'improved' if lo > 0 else 'worse' if hi < 0 else 'uncertain'
    return dict(delta=float(d.mean()), ci95=[lo, hi], receptors=int(len(d)), label=label)


def slice_decisions(wide, info, setting):
    st = SETTINGS[setting]
    core = [n for n, c in info.items() if c['group'] in ('core', 'perm')]
    failed = [n for n in core if info[n]['error']]
    common = wide[[n for n in core if n in wide]].dropna()
    M = common.mean()
    out = dict(receptors_common=int(len(common)), failed_cases=failed, M={k: float(M[k]) for k in CORE if k in M},
               undetermined=any(info[n]['error'] for n in CORE if n in info))   # plan 4: failed core case
    for s in ('raw', 'snorm'):
        perms = sorted(k for k in M.index if k.startswith(f'{s}_perm_'))
        if perms and s in M:
            mp = M[perms].to_numpy()
            p = float((1 + (mp >= M[s] - 1e-12).sum()) / (1 + len(mp)))
            out[f'D1_{s}'] = dict(M=float(M[s]), n_perm=len(mp), perm_mean=float(mp.mean()), perm_sd=float(mp.std(ddof=1)) if len(mp) > 1 else None,
                                  p=p, label='information' if p <= ALPHA else 'no_information')
        out[f'D2_{s}'] = paired(common, s, 'no_ppi')
    out['D2_uniform'] = paired(common, 'uniform', 'no_ppi')
    out['D3_raw_vs_snorm'] = paired(common, 'raw', 'snorm')
    out['D3_raw_vs_snorm']['label'] = {'improved': 'raw_better', 'worse': 'snorm_better', 'uncertain': 'indistinguishable'}.get(
        out['D3_raw_vs_snorm']['label'], out['D3_raw_vs_snorm']['label'])
    sens = {}
    for n, c in info.items():
        if c['group'] == 'sensitivity':
            ref = 'no_ppi_scored_only' if c['panel'] == 'scored' and n != 'no_ppi_scored_only' else 'no_ppi'
            if n != 'no_ppi_scored_only':
                sens[n] = dict(paired(wide, n, ref), reference=ref)
    out['sensitivity'] = sens
    return out


def conclusion(d1, d2):
    if 'undetermined' in (d1, d2):
        return '无法判定(核心条件求解失败)'
    if d1 is None:
        return {'improved': '改进(无 D1)', 'worse': '变差(无 D1)', 'uncertain': '不确定(无 D1)'}.get(d2, d2)
    if d1 == 'information':
        return '在本评价下改进了推断' if d2 == 'improved' else '带有配对信息,但目前不足以改进推断'
    return '在本评价下看不出配对信息'


def stability(folder, pairs=(('raw', 'raw_n200k'), ('snorm', 'snorm_n200k'))):
    """stage2_plan.md 8.7 criteria: per receptor relative L1 <= 1e-3 and Spearman >= 0.999; total transport 1e-3."""
    from scipy.stats import spearmanr
    out = {}
    for short, long_ in pairs:
        A, B = Path(folder)/short, Path(folder)/long_
        if not (A/'signal.npz').exists() or not (B/'signal.npz').exists():
            out[short] = dict(passed=False, note='missing'); continue
        za, zb = np.load(A/'signal.npz'), np.load(B/'signal.npz')
        check(list(za['receptors']) == list(zb['receptors']), 'Receptor order differs')
        rel, rho, zero = [], [], 0
        for x, y in zip(za['signal'], zb['signal']):
            if not y.any() and not x.any():
                zero += 1; continue
            rel.append(np.abs(x - y).sum() / max(np.abs(y).sum(), 1e-300))
            rho.append(spearmanr(x, y)[0] if np.ptp(x) > 0 and np.ptp(y) > 0 else np.nan)
        ta = json.loads((A/'case.json').read_text())['total_transport']; tb = json.loads((B/'case.json').read_text())['total_transport']
        tt = abs(ta - tb) / tb
        passed = bool(max(rel) <= 1e-3 and np.nanmin(rho) >= 0.999 and not np.isnan(rho).any() and tt <= 1e-3)
        out[short] = dict(receptors=len(rel), zero_receptors=zero, max_rel_l1=float(max(rel)), min_spearman=float(np.nanmin(rho)),
                          undefined_spearman=int(np.isnan(rho).sum()), total_transport_rel=float(tt), passed=passed)
    return out


def summarise(args):
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    out.mkdir(parents=True)
    result = {}
    for setting in ('main', 'second'):
        res = {}
        for slice_ in ('post', 'ant'):
            folder = Path(args.runs)/f'{slice_}_{setting}'
            if not (folder/'cases.json').exists():
                continue
            wide, info = load_runs(folder)
            wide.to_csv(out/f'{slice_}_{setting}_receptor_partial_spearman.csv')
            res[slice_] = slice_decisions(wide, info, setting)
            for k in TOP_TF[1:]:                                   # sensitivity: other TF counts
                w2, _ = load_runs(folder, k)
                res[slice_][f'top_tf_{k}'] = {s: paired(w2[[c for c in CORE if c in w2]].dropna(), s, 'no_ppi') for s in ('raw', 'snorm')}
            if setting == 'second' and slice_ == 'post':
                res['stability'] = stability(folder)
        if not res:
            continue
        verdict = {}
        for s in ('raw', 'snorm'):
            post, ant = res.get('post', {}), res.get('ant', {})
            d1p, d1a = (post.get(f'D1_{s}') or {}).get('label'), (ant.get(f'D1_{s}') or {}).get('label')
            d2p, d2a = post.get(f'D2_{s}', {}).get('label'), ant.get(f'D2_{s}', {}).get('label')
            if post.get('undetermined'):
                d1p = d2p = 'undetermined'
            if ant.get('undetermined'):
                d1a = d2a = 'undetermined'
            verdict[s] = dict(D1_post=d1p, D1_ant=d1a, D1_confirmed=bool(d1p is not None and d1p == d1a),
                              D2_post=d2p, D2_ant=d2a, D2_confirmed=bool(d2p is not None and d2p == d2a),
                              conclusion_post=conclusion(d1p, d2p),
                              conclusion_ant=conclusion(d1a, d2a) if ant else None)
        d3p, d3a = res.get('post', {}).get('D3_raw_vs_snorm', {}).get('label'), res.get('ant', {}).get('D3_raw_vs_snorm', {}).get('label')
        verdict['D3'] = dict(post=d3p, ant=d3a, confirmed=bool(d3p is not None and d3p == d3a))
        if setting == 'second':
            stab = res.get('stability', {})
            verdict['stability_passed'] = bool(stab) and all(v.get('passed') for v in stab.values())
            if not verdict['stability_passed']:
                verdict['note'] = 'stability diagnostic failed or missing: second setting undetermined (plan section 2)'
        failed = [f for v in res.values() if isinstance(v, dict) for f in v.get('failed_cases', [])]
        verdict['failed_cases'] = failed
        result[setting] = dict(details=res, verdict=verdict)
    write_json(out/'decisions.json', dict(created_utc=now(), n_boot=N_BOOT, alpha=ALPHA, **result))
    print(json.dumps({k: v['verdict'] for k, v in result.items()}, indent=2, ensure_ascii=False))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    c = sub.add_parser('cases'); c.add_argument('--setting', choices=list(SETTINGS), required=True)
    c.add_argument('--slice', choices=['post', 'ant'], required=True); c.add_argument('--output', required=True)
    c.add_argument('--smoke', action='store_true')
    r = sub.add_parser('run')
    for k in ('input', 'designs', 'ppi-scores', 'protein-map', 'genes', 'db', 'cases', 'output'):
        r.add_argument('--' + k, required=True)
    r.add_argument('--setting', choices=list(SETTINGS), required=True); r.add_argument('--slice', choices=['post', 'ant'], required=True)
    r.add_argument('--case', nargs='+', required=True); r.add_argument('--device', help='cpu | cuda (default: per setting)')
    r.add_argument('--save-case', help='also write stage-2 result folders here (equivalence gate)')
    e = sub.add_parser('equivalence'); e.add_argument('--runs', required=True); e.add_argument('--eval-dir', required=True)
    e.add_argument('--case', nargs='+', required=True); e.add_argument('--output', required=True)
    s = sub.add_parser('summarise'); s.add_argument('--runs', required=True); s.add_argument('--output', required=True)
    return p


if __name__ == '__main__':
    a = parser().parse_args()
    dict(cases=cases_cmd, run=run, equivalence=equivalence, summarise=summarise)[a.command](a)
