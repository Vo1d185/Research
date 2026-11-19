"""Semi-synthetic injection test: can the OT cost make use of a PPI score, and how accurate must it be?
Plan: reports/stage2_inject_01/plan.md. Environment: commot-research. No PDE: the ground truth is an
explicit formula on the real slice (real spot coordinates and real ligand / receptor expression).

  truth      prepared stage-2 input -> known binding constants K (some pairs are decoys, K = 0) and the
             true complex C*_p(y) of every pair p at every spot y:
               ligand reaching y:  L_i(y) = sum_x G(d(x,y)) a_i(x),  G = exp(-d/lambda) (diffusion, within the
                                   cutoff, d = 0 included) or 1 on the contact graph (self + hex neighbours)
               competitive binding: C*_ij(y) = b_j(y) K_ij L_i(y) / (1 + sum_i' K_i'j L_i'(y))
             G is deliberately not the OT kernel F = r0/(r0+d).
  scores     truth -> score vectors in (0,1): latent u (true pairs: 1 + 0.5 z(log K), z clipped to +-2;
             decoys: -1), S = sigmoid(gamma (u + sigma eta)); sigma set so that AUROC(true vs decoy)
             equals each pre-registered level; every vector also gets a within-type permutation
  run        one (score variant, normalisation) over the beta grid:
             cost = -log F(d) - beta * log S_eff, S_eff = S (raw) or S / sum over the ligand's receptors;
             metrics against the truth; transport plans are not saved
  summarise  all runs -> tables and the pre-registered decisions (plan section 5)

Every path is a command-line argument; existing outputs are never overwritten.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.special import expit, ndtr
from scipy.stats import rankdata

GAMMA = 4.0
AUROC_LEVELS = {'q100': 1.0, 'q90': 0.9, 'q80': 0.8, 'q70': 0.7, 'q60': 0.6}
BETAS = (0.05, 0.1, 0.2, 0.5, 1.0)
CONVERGENCE_RULES = ('cap', 'checks')       # plan.md / amendment_02.md (Stage 2 plan section 8.7)
S_FLOOR = 1e-12


def check(ok, msg):
    if not ok:
        raise ValueError(msg)


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as s:
        for block in iter(lambda: s.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def clean(obj):
    """NaN / inf -> None, recursively (JSON has no NaN)."""
    if isinstance(obj, dict):
        return {k: clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, (float, np.floating)):
        return float(obj) if np.isfinite(obj) else None
    return obj


def write_json(path, obj):
    Path(path).write_text(json.dumps(clean(obj), indent=2, ensure_ascii=False, allow_nan=False, default=str))


def fresh(path):
    path = Path(path)
    check(not path.exists(), f'Output exists: {path}')
    path.mkdir(parents=True)
    return path


def load(source):
    import stage2
    prep = stage2.verified(source)
    cfg, x, spatial, lr, ligands, receptors, a, b = stage2.load_input(source)
    edges = dict(zip(['diffusion', 'contact'], stage2.spatial_edges(spatial, cfg['diffusion_cutoff_um'])))
    return prep, cfg, x, spatial, lr.reset_index(drop=True), ligands, receptors, a, b, edges


# ------------------------------------------------------------------ truth
def ligand_fields(x, lr, edges, lam):
    """L_i(y) for every ligand, using the edges of its (single) interaction type."""
    n = len(x)
    types = lr.groupby('ligand').type.unique()
    check((types.map(len) == 1).all(), 'A ligand has pairs of both types')
    out = {}
    for lig, t in types.items():
        s, tt, d = edges[t[0]]
        g = np.exp(-d / lam) if t[0] == 'diffusion' else np.ones(len(d))
        out[lig] = np.bincount(tt, weights=g * x[lig].to_numpy(float)[s], minlength=n)
    return out


def complexes(lr, K, L, x):
    """C*_p(y) for every pair (rows follow lr), competitive binding at each receptor."""
    n = len(x)
    load_ = {r: np.zeros(n) for r in lr.receptor.unique()}
    for p, r in enumerate(lr.itertuples(index=False)):
        load_[r.receptor] += K[p] * L[r.ligand]
    C = np.zeros((len(lr), n))
    for p, r in enumerate(lr.itertuples(index=False)):
        C[p] = x[r.receptor].to_numpy(float) * K[p] * L[r.ligand] / (1.0 + load_[r.receptor])
    return C, load_


def truth(args):
    prep, cfg, x, spatial, lr, ligands, receptors, a, b, edges = load(args.input)
    rng = np.random.default_rng(args.seed)
    decoy = np.zeros(len(lr), bool)
    for t, idx in lr.groupby('type').groups.items():
        idx = np.asarray(sorted(idx))
        decoy[rng.choice(idx, size=int(round(args.decoy_fraction * len(idx))), replace=False)] = True
    u = rng.uniform(-args.log10k_range / 2, args.log10k_range / 2, len(lr))
    L = ligand_fields(x, lr, edges, args.lambda_um)
    rel = np.where(decoy, 0.0, 10.0 ** u)
    _, z = complexes(lr, rel, L, x)
    expressed = np.concatenate([z[r][x[r].to_numpy(float) > 0] for r in z])
    check((expressed > 0).any(), 'No receptor spot receives ligand')
    kappa = 1.0 / float(np.median(expressed[expressed > 0]))
    K = rel * kappa
    C, load_ = complexes(lr, K, L, x)
    out = fresh(args.output)
    pairs = lr[['ligand', 'receptor', 'type']].assign(decoy=decoy, log10K_rel=np.where(decoy, np.nan, u), K=K,
                                                     C_total=C.sum(axis=1))
    pairs.to_csv(out/'pairs.csv', index=False)
    np.savez_compressed(out/'truth.npz', C=C, K=K, decoy=decoy)
    occ = np.concatenate([load_[r][x[r].to_numpy(float) > 0] for r in load_])
    write_json(out/'truth.json', dict(created_utc=now(), input=str(Path(args.input).resolve()),
        input_sha256=prep['prepared_sha256'], seed=args.seed, lambda_um=args.lambda_um,
        decoy_fraction=args.decoy_fraction, log10k_range=args.log10k_range, kappa=kappa,
        n_pairs=len(lr), n_decoys=int(decoy.sum()), n_spots=len(x),
        occupancy_at_receptor_spots=dict(median=float(np.median(occ)), q90=float(np.quantile(occ, .9)),
                                         frac_gt_1=float((occ > 1).mean())),
        truth_sha256=digest(out/'truth.npz'),
        model='L_i(y)=sum_x G(d) a_i(x); C*=b_j K L_i/(1+sum K L); ligand depletion ignored'))
    print(f'truth: {len(lr)} pairs ({int(decoy.sum())} decoys), kappa={kappa:.3g} -> {out}')


# ------------------------------------------------------------------ scores
def latent(K, decoy):
    t = np.log10(np.where(decoy, 1.0, K))
    zt = (t[~decoy] - t[~decoy].mean()) / (t[~decoy].std() or 1.0)
    u = np.full(len(K), -1.0)
    u[~decoy] = 1.0 + 0.5 * np.clip(zt, -2, 2)
    return u


def expected_auroc(u, decoy, sigma):
    """P(u_T + sigma e1 > u_D + sigma e2) averaged over true x decoy pairs (closed form)."""
    diff = u[~decoy][:, None] - u[decoy][None, :]
    return float(ndtr(diff / (sigma * np.sqrt(2))).mean()) if sigma > 0 else float((diff > 0).mean() + 0.5 * (diff == 0).mean())


def sigma_for(u, decoy, target):
    if target >= 1.0:
        return 0.0
    lo, hi = 1e-6, 100.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if expected_auroc(u, decoy, mid) > target else (lo, mid)
    return 0.5 * (lo + hi)


def auroc(pos, neg):
    r = rankdata(np.concatenate([pos, neg]))
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def spearman(x, y):
    x, y = rankdata(x), rankdata(y)
    sx, sy = x.std(), y.std()
    return float(((x - x.mean()) * (y - y.mean())).mean() / (sx * sy)) if sx > 0 and sy > 0 else np.nan


def scores(args):
    tdir = Path(args.truth)
    pairs = pd.read_csv(tdir/'pairs.csv')
    t = np.load(tdir/'truth.npz')
    K, decoy = t['K'], t['decoy']
    rng = np.random.default_rng(args.seed)
    u = latent(K, decoy)
    table = pairs[['ligand', 'receptor', 'type', 'decoy', 'K']].copy()
    info = {}
    for name, level in AUROC_LEVELS.items():
        sigma = sigma_for(u, decoy, level)
        s = expit(GAMMA * (u + sigma * rng.standard_normal(len(u))))
        perm = s.copy()
        for _, idx in pairs.groupby('type').groups.items():
            idx = np.asarray(sorted(idx))
            perm[idx] = s[rng.permutation(idx)]
        table[name], table[name + '_perm'] = s, perm
        info[name] = dict(target_auroc=level, sigma=sigma, achieved_auroc=auroc(s[~decoy], s[decoy]),
                          spearman_with_K_true_pairs=spearman(s[~decoy], K[~decoy]),
                          perm_auroc=auroc(perm[~decoy], perm[decoy]),
                          frac_gt_0_99=float((s > 0.99).mean()), frac_lt_0_01=float((s < 0.01).mean()))
    out = fresh(args.output)
    table.to_csv(out/'scores.csv', index=False)
    write_json(out/'scores.json', dict(created_utc=now(), truth=str(tdir.resolve()),
               truth_sha256=digest(tdir/'truth.npz'), seed=args.seed, gamma=GAMMA, levels=info,
               scores_sha256=digest(out/'scores.csv')))
    print(json.dumps({k: {kk: round(v[kk], 3) for kk in ('sigma', 'achieved_auroc', 'spearman_with_K_true_pairs')}
                      for k, v in info.items()}, indent=1))


# ------------------------------------------------------------------ run
def effective_score(lr, s, norm):
    s = np.clip(np.asarray(s, float), S_FLOOR, 1.0)
    if norm == 'raw':
        return s
    check(norm == 'ligand', f'Unknown normalisation {norm}')
    tot = pd.Series(s).groupby(lr.ligand.to_numpy()).transform('sum').to_numpy()
    return s / tot


def received(T, col, p, n, n_pairs):
    """r_hat_p(y) = sum_x T_p(x, y): column index = receptor * n + spot."""
    out = np.zeros((n_pairs, n))
    np.add.at(out, (p, col % n), T)
    return out


def metrics(R, C, lr, x, decoy, T, p):
    expressed = np.stack([x[r].to_numpy(float) > 0 for r in lr.receptor])          # pair x spot
    pooled = spearman(R[expressed], C[expressed])
    rec = []
    for r, idx in lr.groupby('receptor').groups.items():
        idx = np.asarray(sorted(idx)); m = x[r].to_numpy(float) > 0
        rec.append(spearman(R[idx].sum(0)[m], C[idx].sum(0)[m]))
    rec = np.asarray(rec, float)
    tv = []
    for _, idx in lr[~decoy].groupby('ligand').groups.items():
        idx = np.asarray(sorted(idx))
        if len(idx) < 2:
            continue
        w_hat, w = R[idx].sum(1), C[idx].sum(1)
        if w_hat.sum() > 0 and w.sum() > 0:
            tv.append(0.5 * np.abs(w_hat / w_hat.sum() - w / w.sum()).sum())
    total = float(T.sum())
    return dict(pooled_spearman=pooled, receptor_spearman_mean=float(np.nanmean(rec)),
                receptor_spearman_n=int(np.isfinite(rec).sum()),
                ligand_allocation_tv_mean=float(np.mean(tv)) if tv else np.nan, ligand_allocation_n=len(tv),
                decoy_transport_fraction=float(T[decoy[p]].sum() / total) if total else np.nan,
                total_transport=total)


def run(args):
    import stage2
    from research import solve_sparse
    device = getattr(args, 'device', 'cpu') or 'cpu'
    if device != 'cpu':                       # amendment_01.md: torch port of the same solver, same checks
        from functools import partial
        from gpu_unot import solve_sparse_gpu
        solve_sparse = partial(solve_sparse_gpu, device=device)
    prep, cfg, x, spatial, lr, ligands, receptors, a, b, edges = load(args.input)
    cfg = stage2.run_config(cfg, args)
    tdir, sdir = Path(args.truth), Path(args.scores)
    t = np.load(tdir/'truth.npz')
    C, decoy = t['C'], t['decoy']
    pairs = pd.read_csv(tdir/'pairs.csv')
    check(pairs[['ligand', 'receptor']].equals(lr[['ligand', 'receptor']]), 'Truth was built from another input')
    sc = pd.read_csv(sdir/'scores.csv')
    check(sc[['ligand', 'receptor']].equals(lr[['ligand', 'receptor']]), 'Scores do not match the input')
    n = len(x)
    row, col, c0, d, p = stage2.build_sparse_cost(lr, ligands, receptors, a, b, n, edges, args.r0_um)
    if args.variant == 'none':
        betas, s_eff = [0.0], np.ones(len(lr))
    elif args.variant == 'uniform':
        betas = list(args.beta)
        s_eff = 1.0 / lr.groupby('ligand').receptor.transform('size').to_numpy(float)
    else:
        check(args.variant in sc, f'Unknown score variant {args.variant}')
        betas, s_eff = list(args.beta), effective_score(lr, sc[args.variant], args.norm)
    out = fresh(args.output)
    rows = []
    for beta in betas:
        c = c0 + beta * (-np.log(s_eff))[p]
        start = time.time()
        try:
            T, info = solve_sparse(row, col, c, a, b, cfg['eps'], cfg['rho'], cfg['nitermax'], cfg['stopthr'])
            error = ''
        except ValueError as e:                     # mass check failed: not converged within nitermax
            check('not converged' in str(e), str(e))
            T, info, error = None, dict(cap_reached=True), str(e)
        seconds = time.time() - start
        if T is None:
            m = {k: np.nan for k in ('pooled_spearman', 'receptor_spearman_mean', 'ligand_allocation_tv_mean',
                                     'decoy_transport_fraction', 'total_transport')}
            m.update(receptor_spearman_n=0, ligand_allocation_n=0)
        else:
            m = metrics(received(T, col, p, n, len(lr)), C, lr, x, decoy, T, p)
        rows.append(dict(seed=args.seed, device=device, variant=args.variant, norm=args.norm if args.variant not in ('none', 'uniform') else args.variant,
                         beta=beta, eps=cfg['eps'], rho=cfg['rho'], r0_um=args.r0_um, niter=info.get('niter'),
                         cap_reached=info.get('cap_reached'), solve_seconds=round(seconds, 1),
                         frac_ligand=float(T.sum() / a.sum()) if T is not None else np.nan, error=error, **m,
                         **{f'check_{k}': v for k, v in info.items() if k.startswith(('mass', 'kkt', 'max'))}))
        print(f'{args.variant}/{args.norm} beta={beta}: pooled={m["pooled_spearman"]:.4f} '
              f'decoy={m["decoy_transport_fraction"]:.3f} {seconds:.0f}s {info.get("niter")} it'
              + (' NOT CONVERGED' if error else ''), flush=True)
    pd.DataFrame(rows).to_csv(out/'metrics.csv', index=False)
    write_json(out/'run.json', dict(created_utc=now(), input=str(Path(args.input).resolve()),
        truth_sha256=digest(tdir/'truth.npz'), scores_sha256=digest(sdir/'scores.csv'), config=cfg,
        variant=args.variant, norm=args.norm, betas=betas, r0_um=args.r0_um, entries=len(c0), device=device))


# ------------------------------------------------------------------ summarise
def converged(m, rule='cap'):
    """cap (plan.md): converged = iteration cap not reached.
    checks (amendment_02.md = Stage 2 plan section 8.7): converged = the solver's mass and KKT checks passed (the run
    has a transport plan and no error); reaching the cap is recorded but is not a failure."""
    check(rule in CONVERGENCE_RULES, f'Unknown convergence rule {rule}')
    if rule == 'cap':
        return ~m.cap_reached.fillna(False).astype(bool)
    err = m['error'].fillna('').astype(str).str.len() > 0 if 'error' in m else pd.Series(False, index=m.index)
    return ~err & m.total_transport.notna()


def decide(m, v12=None, allow_unconverged=False, rule='cap'):
    """Pre-registered decisions (plan section 5) from one parameter setting's metrics table."""
    key = ['seed', 'variant', 'norm', 'beta']
    check(not m.duplicated(key).any(), 'Duplicate runs')
    m = m.copy()
    m['converged'] = converged(m, rule)
    base = m[m.variant == 'none'].set_index('seed')
    check(allow_unconverged or base.converged.all(), 'A no-PPI baseline did not converge')
    m['gain'] = m.pooled_spearman - m.seed.map(base.pooled_spearman)
    m['gain_receptor'] = m.receptor_spearman_mean - m.seed.map(base.receptor_spearman_mean)
    real = m[m.variant.isin(AUROC_LEVELS)].copy()
    perm = m[m.variant.str.endswith('_perm')].copy()
    perm['variant'] = perm.variant.str.replace('_perm', '', regex=False)
    real = real.merge(perm[key + ['pooled_spearman', 'receptor_spearman_mean', 'converged']], on=key,
                      suffixes=('', '_perm'))
    real['ok'] = real.converged & real.converged_perm
    real['gain_vs_perm'] = real.pooled_spearman - real.pooled_spearman_perm
    real['gain_vs_perm_receptor'] = real.receptor_spearman_mean - real.receptor_spearman_mean_perm
    seeds = sorted(m.seed.unique())
    check(len(seeds) >= 2, 'Need at least two seeds (the first selects, the others confirm)')
    sel, conf = seeds[0], seeds[1:]
    dec = dict(seeds=[int(x) for x in seeds], selection_seed=int(sel), convergence_rule=rule,
               runs=int(len(m)), runs_not_converged=int((~m.converged).sum()),
               runs_at_cap=int(m.cap_reached.fillna(False).astype(bool).sum()))

    def all_seeds(g, cond):
        return bool(g.seed.nunique() == len(seeds) and g.ok.all() and cond(g).all())
    o = real[real.variant == 'q100']
    cells = o.groupby(['norm', 'beta'])
    r1 = {f'{k[0]}|{k[1]:g}': all_seeds(g, lambda x: (x.gain > 0) & (x.gain_vs_perm > 0)) for k, g in cells}
    dec['R1_oracle_usable'] = any(r1.values()); dec['R1_cells'] = r1
    eps = float(m.eps.iloc[0])
    raw = o[(o.norm == 'raw') & o.ok].groupby('beta').pooled_spearman.mean()
    dec['R2_best_beta_oracle_raw'] = float(raw.idxmax()) if len(raw) else None
    dec['R2_near_eps'] = bool(len(raw) and np.isclose(raw.idxmax(), [eps / 2, eps, 2 * eps]).any())
    thr = {}
    for v, level in AUROC_LEVELS.items():
        g = real[(real.variant == v) & real.ok]
        s1 = g[g.seed == sel]
        if not len(s1):
            thr[v] = dict(auroc=level, usable=False, note='no converged run on the selection seed'); continue
        pick = s1.sort_values(['gain_vs_perm', 'beta'], ascending=[False, True]).iloc[0]
        c = g[g.seed.isin(conf) & (g.norm == pick.norm) & np.isclose(g.beta, pick.beta)]
        thr[v] = dict(auroc=level, selected_norm=pick.norm, selected_beta=float(pick.beta),
                      selection_gain_vs_perm=float(pick.gain_vs_perm),
                      confirm_gain_vs_perm=[float(x) for x in c.gain_vs_perm],
                      confirm_gain_vs_none=[float(x) for x in c.gain],
                      usable=bool(len(c) == len(conf) and (c.gain_vs_perm > 0).all()))
    dec['R3_by_level'] = thr
    usable = [AUROC_LEVELS[v] for v, r in thr.items() if r['usable']]
    dec['R3_threshold_auroc'] = min(usable) if usable else None
    order = sorted(AUROC_LEVELS, key=AUROC_LEVELS.get)
    flags = [thr[v]['usable'] for v in order]
    dec['R3_monotone'] = bool(all(a <= b for a, b in zip(flags, flags[1:])))
    r5 = {f'{k[0]}|{k[1]:g}': all_seeds(g, lambda x: (x.gain_receptor > 0) & (x.gain_vs_perm_receptor > 0))
          for k, g in cells}
    dec['R5_receptor_metric_sees_oracle'] = any(r5.values()); dec['R5_cells'] = r5
    if v12 is not None:
        t = dec['R3_threshold_auroc']
        lo, hi = v12['ci95']
        dec['R4_v12'] = dict(auroc=v12['auroc'], ci95=[lo, hi], threshold=t,
            placement=('no threshold (no level usable)' if t is None else
                       'above' if lo >= t else 'below' if hi < t else 'undetermined'))
    u = m[m.variant == 'uniform'].groupby('beta')[['gain', 'gain_receptor']].mean()
    dec['uniform_mean_gain_by_beta'] = {f'{k:g}': dict(pooled=float(r.gain), receptor=float(r.gain_receptor))
                                        for k, r in u.iterrows()}
    return dec, real


def single_device(m):
    """One summary = one solver backend (amendment_01.md): never pool CPU and GPU runs."""
    devices = sorted(set(m['device'].fillna('cpu'))) if 'device' in m else ['cpu']
    check(len(devices) == 1, f'Runs from more than one device: {devices}')
    return devices[0]


def summarise(args):
    frames = [pd.read_csv(f) for f in sorted(Path(args.runs).glob('**/metrics.csv'))]
    check(frames, 'No metrics.csv found')
    m = pd.concat(frames, ignore_index=True)
    device = single_device(m)
    v12 = None
    if args.v12_report:
        rep = json.loads(Path(args.v12_report).read_text())
        v12 = (rep.get('panel_vs_unannotated') or {}).get('all')
        check(v12 and v12.get('ci95'), 'V12 report has no panel_vs_unannotated AUROC')
    out = fresh(args.output)
    m.to_csv(out/'all_metrics.csv', index=False)
    dec, real = decide(m, v12, args.allow_unconverged, args.convergence)
    real.to_csv(out/'real_vs_perm.csv', index=False)
    real.groupby(['variant', 'norm', 'beta'])[['gain', 'gain_vs_perm', 'gain_receptor', 'gain_vs_perm_receptor',
                                              'decoy_transport_fraction']].mean().round(5).to_csv(out/'mean_table.csv')
    write_json(out/'decisions.json', dict(created_utc=now(), eps=float(m.eps.iloc[0]), rho=float(m.rho.iloc[0]),
                                          device=device, v12_report=args.v12_report, **dec))
    print(json.dumps({k: dec.get(k) for k in ['R1_oracle_usable', 'R2_best_beta_oracle_raw', 'R2_near_eps',
                                              'R3_threshold_auroc', 'R3_monotone', 'R4_v12',
                                              'R5_receptor_metric_sees_oracle', 'runs_not_converged', 'runs_at_cap',
                                              'convergence_rule']}, indent=2))


STABILITY_TOL = 1e-3          # Stage 2 plan section 8.7: relative L1 / total transport 1e-3
STABILITY_METRICS = ['pooled_spearman', 'receptor_spearman_mean', 'ligand_allocation_tv_mean',
                     'decoy_transport_fraction', 'frac_ligand']
# amendment_02.md section 3: representative main-setting runs repeated with a longer cap (seed 1 only: the seed changes
# the synthetic scores, not the OT budgets). Covers none, uniform, oracle and its permutation, both normalisations,
# the best and the worst score quality, every beta.
STABILITY_RUNS = ('s1_none', 's1_uniform', 's1_q100_raw', 's1_q100_ligand', 's1_q100_perm_raw', 's1_q100_perm_ligand',
                  's1_q60_raw', 's1_q60_perm_ligand')


def run_metrics(folder):
    frames = [pd.read_csv(f) for f in sorted(Path(folder).glob('*/metrics.csv'))]
    check(frames, f'No metrics.csv under {folder}')
    return pd.concat(frames, ignore_index=True)


def stability_check(used_m, long_m):
    """amendment_02.md (Stage 2 section 8.7): every run repeated with the longer cap must (1) pass the solver checks at
    both caps, (2) give metrics within STABILITY_TOL of the used result (total transport: relative). Differences only."""
    key = ['seed', 'variant', 'norm', 'beta']
    check(not long_m.duplicated(key).any() and not used_m.duplicated(key).any(), 'Duplicate runs')
    Lg = long_m.set_index(key).sort_index()
    U = used_m.set_index(key)
    check(Lg.index.isin(U.index).all(), 'A longer-cap run has no counterpart among the used runs')
    U = U.loc[Lg.index]
    failed_used = int((~converged(U.reset_index(), 'checks')).sum())
    failed_long = int((~converged(Lg.reset_index(), 'checks')).sum())
    diffs, nan_mismatch = {}, {}
    for c in STABILITY_METRICS:
        nan_mismatch[c] = int((U[c].isna() != Lg[c].isna()).sum())
        d = (U[c] - Lg[c]).abs()
        diffs[c] = float(d.max()) if d.notna().any() else None
    tt = (U.total_transport - Lg.total_transport).abs() / Lg.total_transport.abs()
    diffs['total_transport_relative'] = float(tt.max()) if tt.notna().any() else None
    within = all(v is not None and v <= STABILITY_TOL for v in diffs.values())
    passed = bool(failed_used == 0 and failed_long == 0 and within and not any(nan_mismatch.values()))
    return dict(compared_solves=int(len(Lg)), tolerance=STABILITY_TOL, failed_checks_used=failed_used,
                failed_checks_long=failed_long, max_abs_diff=diffs, nan_mismatch=nan_mismatch,
                niter_long_max=float(Lg.niter.max()) if 'niter' in Lg and Lg.niter.notna().any() else None,
                passed=passed,
                verdict='accept all main-setting results' if passed else 'main setting undetermined (amendment_02 section 3)')


def stability(args):
    res = stability_check(run_metrics(args.used), run_metrics(args.long))
    res.update(created_utc=now(), used=str(args.used), long=str(args.long))
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    write_json(out, res)
    print(json.dumps({k: res[k] for k in ('compared_solves', 'failed_checks_used', 'failed_checks_long', 'max_abs_diff',
                                          'passed', 'verdict')}, indent=2))


def recheck_commands(args):
    """amendment_02.md: selected main-setting lines of an existing commands.txt, unchanged except a longer --nitermax and
    outputs/logs under --check-runs (truth and scores are reused read-only). Refuses if any original run is missing."""
    import re
    runs, chk = str(args.runs).rstrip('/'), str(args.check_runs).rstrip('/')
    lines = [l for l in Path(args.commands).read_text().splitlines() if l.strip()]
    missing, out, niters, devices = [], [], set(), set()
    for line in lines:
        o = re.search(r'--output (\S+)', line)
        check(o is not None, f'No --output in: {line[:80]}')
        if not (Path(o.group(1))/'metrics.csv').exists():
            missing.append(o.group(1))
        if f'--output {runs}/main/' not in line or Path(o.group(1)).name not in args.select:
            continue
        check('--eps 0.1 --rho 10' in line, f'Main-setting line without --eps 0.1 --rho 10: {line[:80]}')
        n = re.search(r'--nitermax (\d+)', line)
        check(n is not None, 'No --nitermax in a main-setting line')
        niters.add(int(n.group(1)))
        d = re.search(r'--device (\S+)', line); devices.add(d.group(1) if d else 'cpu')
        new = line
        for a, b in ((f'--nitermax {n.group(1)}', f'--nitermax {args.nitermax}'),
                     (f'--output {runs}/main/', f'--output {chk}/main/'), (f'> {runs}/logs/', f'> {chk}/logs/')):
            check(new.count(a) == 1, f'Expected exactly one "{a}" in: {line[:80]}')
            new = new.replace(a, b)
        out.append(new)
    check(not missing, f'{len(missing)} original runs have no metrics.csv (first: {missing[:1]})')
    check(out, f'No main-setting lines for --runs {runs}; give it exactly as in the original run')
    check(len(niters) == 1 and args.nitermax > max(niters), f'Check cap must exceed the original cap {sorted(niters)}')
    found = sorted(Path(re.search(r'--output (\S+)', l).group(1)).name for l in out)
    check(not args.require_all or found == sorted(args.select), f'Missing stability runs: {sorted(set(args.select) - set(found))}')
    check(len(devices) == 1, f'Main-setting runs used more than one device: {sorted(devices)}')
    dst = Path(args.output)
    check(not dst.exists(), f'{dst} exists')
    dst.write_text('\n'.join(out) + '\n')
    info = dict(lines_total=len(lines), selected_lines=len(out), selected=found, original_nitermax=min(niters),
                check_nitermax=args.nitermax, device=devices.pop())
    print(json.dumps(info))
    return info


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    t = sub.add_parser('truth'); t.add_argument('--input', required=True); t.add_argument('--output', required=True)
    t.add_argument('--seed', type=int, required=True)
    t.add_argument('--lambda-um', type=float, default=100.0)
    t.add_argument('--decoy-fraction', type=float, default=1 / 3)
    t.add_argument('--log10k-range', type=float, default=2.0)
    s = sub.add_parser('scores'); s.add_argument('--truth', required=True); s.add_argument('--output', required=True)
    s.add_argument('--seed', type=int, required=True)
    r = sub.add_parser('run'); r.add_argument('--input', required=True); r.add_argument('--truth', required=True)
    r.add_argument('--scores', required=True); r.add_argument('--output', required=True)
    r.add_argument('--seed', type=int, required=True)
    r.add_argument('--variant', required=True, help="none | uniform | q100 | q100_perm | q90 | ...")
    r.add_argument('--norm', choices=['raw', 'ligand'], default='raw')
    r.add_argument('--beta', type=float, nargs='+', default=list(BETAS))
    r.add_argument('--r0-um', type=float, default=50.0)
    for k in ('eps', 'rho', 'diffusion_cutoff_um', 'stopthr'):
        r.add_argument('--' + k.replace('_', '-'), type=float)
    r.add_argument('--nitermax', type=int)
    r.add_argument('--device', default='cpu', help='cpu (official solver) | cuda, cuda:N (gpu_unot.py port, amendment_01) | torch-cpu (port on CPU, tests only)')
    m = sub.add_parser('summarise'); m.add_argument('--runs', required=True); m.add_argument('--output', required=True)
    m.add_argument('--v12-report', help='reports/stage2_ppi_01/scores/score_report.json (R4)')
    m.add_argument('--allow-unconverged', action='store_true', help='smoke test only')
    m.add_argument('--convergence', choices=CONVERGENCE_RULES, default='cap',
                   help='cap: plan.md (cap reached = not converged); checks: amendment_02.md (Stage 2 section 8.7)')
    st = sub.add_parser('stability', help='amendment_02.md: used runs vs the same runs with a longer cap')
    st.add_argument('--used', required=True, help='<runs>/main'); st.add_argument('--long', required=True, help='<check-runs>/main')
    st.add_argument('--output', required=True)
    rc = sub.add_parser('recheck-commands', help='amendment_02.md: main-setting commands again at a lower cap')
    for k in ('commands', 'runs', 'check-runs', 'output'):
        rc.add_argument('--' + k, required=True)
    rc.add_argument('--nitermax', type=int, required=True)
    rc.add_argument('--select', nargs='+', default=list(STABILITY_RUNS))
    rc.add_argument('--require-all', action='store_true', help='main run: every STABILITY_RUNS entry must be present')
    return p


if __name__ == '__main__':
    a = parser().parse_args()
    dict(truth=truth, scores=scores, run=run, summarise=summarise, stability=stability,
         **{'recheck-commands': recheck_commands})[a.command](a)
