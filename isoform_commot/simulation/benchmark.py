"""PDE r0-lambda benchmark (reports/pde_r0_plan.md).

  python -m simulation.benchmark calibrate --archive ZIP --output reports/pde_r0_calibration
  python -m simulation.benchmark scenarios --archive ZIP --calibration reports/pde_r0_calibration --output DIR
  python -m simulation.benchmark truth     --scenarios DIR --scenario NAME --output DIR [--dynamic]
  python -m simulation.benchmark ot        --scenarios DIR --scenario NAME --k K --condition COND --eps E --rho R
                                           [--nitermax N] --output DIR
  python -m simulation.benchmark jobs      --scenarios DIR --results DIR --output JOBS.tsv
  python -m simulation.benchmark pipeline  --archive ZIP --results DIR --output DIR
  python -m simulation.benchmark stability --results DIR --output DIR
  python -m simulation.benchmark evaluate  --scenarios DIR --results DIR --output DIR
  python -m simulation.benchmark analyze   --metrics DIR --output DIR

calibrate  archived input statistics and the one-off generator calibration (section 4.2; no PDE, no OT)
scenarios  generated scenarios (structures 2-10 = archive cases, repeats 1-5) and archived scenarios (sim 1)
truth      steady state with checks (5.1 items 1-2), t = 50 (5.2); --dynamic adds the BDF check (5.1 item 3)
ot         one OT run: inv | r0_<r0/lambda> | cd_<absolute dis_thr> | ref (sections 6-7). COMMOT's transport
           depends on the absolute dis_thr only, so cd runs are keyed by it and shared across lambda.
jobs       job list in the priority order of section 11 (for simulation.jobqueue)
pipeline   condition R vs the archive's own COT output (section 10)
stability  D1-D3: 50000 vs 200000 iterations (section 10, stage-2 8.7 thresholds)
evaluate   per-pair Spearman and RMSE of every run against both truths (section 8)
analyze    section 9: r0* selection, hypothesis H, bootstrap, descriptive comparisons (simulation/analysis.py)
"""
import argparse
from datetime import datetime, timezone
from dataclasses import asdict
import hashlib
import io
import json
from pathlib import Path
import time
import zipfile

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import spearmanr

from research import check, launch_info
from simulation import fields
from simulation.pde import Grid, Model, System

ARCHIVE_SHA256 = 'ec8eb74c878f299b0b7c017a1ca758777c188c2ea1c672120e62b2de42ccc50f'
RATES = dict(a=1.0, b=1e-3, c=0.2)
N_GRID, H = 50, 8.0
LAMBDA_A = float(np.sqrt(200.0 / 0.2))
KS = (-2, -1, 0, 1, 2)
R0_OVER_LAMBDA = (0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
DIS_OVER_LAMBDA = (1.0, 2.0, 4.0)
DIS_FIXED = 150.0
CUTOFF_OVER_LAMBDA = 4.0
STRUCTURES = tuple(range(2, 11))
REPEATS = tuple(range(1, 6))
SETTINGS = ((0.1, 10.0), (0.2, 1.0))           # main, sensitivity
NITERMAX, STOPTHR = 50000, 1e-11
REFERENCE = dict(eps=0.25, rho=1.0, dis_thr=150.0, nitermax=10000, stopthr=1e-8)
DIAGNOSTICS = (('gen_c06_s1', -2, 'r0_0.125'), ('gen_c06_s1', 2, 'r0_8'), ('gen_c06_s1', 0, 'cd_4lambda'))
STABILITY = dict(rel_l1=1e-3, spearman=0.999, total=1e-3)
CHECKS = dict(linear=1e-10, nonnegative=1e-12, conservation=1e-12, dynamic_change=1e-8, dynamic_match=1e-5)


def lam(k):
    return LAMBDA_A * 2.0 ** k


def fmt(x):
    return f'{x:.6g}'


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, default=str))


def open_archive(path):
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    check(digest == ARCHIVE_SHA256, f'Unexpected archive checksum {digest}')
    return zipfile.ZipFile(path)


def archived(zf, case, sim):
    from simulation.examples import load_archived
    system, *_ = load_archived(zf, case, sim)
    return system.p, system.r0, (system.m.a > 0).astype(float)


# ------------------------------------------------------------------ calibration and scenarios (4.2, 4.3)
def calibrate(args):
    launch = launch_info()
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    zf = open_archive(args.archive)
    stats, structures = [], {}
    for case in STRUCTURES:
        for sim in range(1, 11):
            p, r, binding = archived(zf, case, sim)
            structures.setdefault(case, binding)
            check(np.array_equal(structures[case], binding), f'case {case}: binding differs between simulations')
            stats.append(fields.scenario_statistics(p, r, binding))
    reference = fields.pooled(stats)
    rows, chosen = fields.calibrate(reference, structures)
    sample = fields.pooled([fields.scenario_statistics(*fields.scenario_fields(c, s, a.shape[0], a.shape[1], chosen), a)
                            for c, a in sorted(structures.items()) for s in fields.CALIBRATION_REPEATS])
    _, parts = fields.ks_objective(reference, sample)
    comparison = []
    for k in fields.STATISTICS:
        for source, v in [('archive', reference[k]), ('generated', sample[k])]:
            comparison.append(dict(statistic=k, source=source, n=len(v), median=float(np.median(v)),
                                   q10=float(np.quantile(v, 0.1)), q90=float(np.quantile(v, 0.9)),
                                   mean=float(v.mean()), ks=parts[k]))
    out.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(out/'calibration_grid.csv', index=False)
    pd.DataFrame(comparison).to_csv(out/'comparison.csv', index=False)
    write_json(out/'calibration.json', dict(status='complete', launch=launch, finished=now(),
        archive=str(args.archive), archive_sha256=ARCHIVE_SHA256, reference='archive cases 2-10, sims 1-10 (inputs only)',
        calibration_seeds='default_rng([structure, s]) for s in 1001..1010', grid=fields.CALIBRATION_GRID,
        statistics=list(fields.STATISTICS), objective='sum of two-sample KS statistics', chosen=asdict(chosen),
        objective_value=sum(parts.values()), ks=parts, structures={c: a.tolist() for c, a in structures.items()}))
    print(json.dumps(dict(chosen=asdict(chosen), ks=parts), indent=2))


def scenarios(args):
    launch = launch_info()
    out = Path(args.output)
    check(not out.exists(), f'Output exists: {out}')
    cal = json.loads((Path(args.calibration)/'calibration.json').read_text())
    params = fields.FieldParams(**cal['chosen'])
    zf = open_archive(args.archive)
    out.mkdir(parents=True)
    names = []
    for case in STRUCTURES:
        p_arc, r_arc, binding = archived(zf, case, 1)
        for s in REPEATS:
            p, r = fields.scenario_fields(case, s, *binding.shape, params)
            name = f'gen_c{case:02d}_s{s}'
            np.savez_compressed(out/f'{name}.npz', p=p, r0=r, binding=binding, structure=case, repeat=s,
                                source='generated')
            names.append(name)
        name = f'arc_c{case:02d}_sim1'
        np.savez_compressed(out/f'{name}.npz', p=p_arc, r0=r_arc, binding=binding, structure=case, repeat=1,
                            source='archive')
        names.append(name)
    write_json(out/'manifest.json', dict(status='complete', launch=launch, finished=now(), params=asdict(params),
        calibration=str(args.calibration), archive_sha256=ARCHIVE_SHA256, scenarios=names,
        seeds='generated: numpy.random.default_rng([structure, repeat]); archive: case c, sim 1 inputs'))


def load_scenario(folder, name):
    with np.load(Path(folder)/f'{name}.npz', allow_pickle=False) as z:
        return {k: z[k] for k in z.files}


def pairs_of(binding):
    return [(i, j) for i in range(binding.shape[0]) for j in range(binding.shape[1]) if binding[i, j] > 0]


# ------------------------------------------------------------------ truth (5.1, 5.2)
def system_for(scen, k):
    binding = scen['binding']
    nl, nr = binding.shape
    D = RATES['c'] * lam(k) ** 2
    model = Model(D=[D] * nl, a=RATES['a'] * binding, b=np.full((nl, nr), RATES['b']), c=[RATES['c']] * nl)
    check(model.pairs == pairs_of(binding), 'Pair order differs')
    return System(model, Grid(N_GRID, -200.0, 200.0, 2), scen['p'], scen['r0'])


def steady_checks(system, L, R, C):
    lap = system.grid.laplacian_matrix()
    linear = []
    for i in range(system.m.nl):
        res = (system.m.c[i] * L[i].ravel() - system.m.D[i] * (lap @ L[i].ravel())) - system.p[i].ravel()
        linear.append(float(np.abs(res).max() / np.abs(system.p[i]).max()))
    nonneg = [float(-min(L[i].min(), 0) / L[i].max()) for i in range(system.m.nl)]
    total = R.copy()
    for k, (i, j) in enumerate(system.m.pairs):
        total[j] = total[j] + C[k]
    conservation = [float(np.abs(total[j] - system.r0[j]).max() / system.r0[j].max()) for j in range(system.m.nr)]
    ok = max(linear) <= CHECKS['linear'] and max(nonneg) <= CHECKS['nonnegative'] and \
        max(conservation) <= CHECKS['conservation']
    return dict(linear_residual=max(linear), negative_ligand=max(nonneg), conservation=max(conservation), passed=ok)


def truth(args):
    launch = launch_info()
    scen = load_scenario(args.scenarios, args.scenario)
    out = Path(args.output)/args.scenario
    out.mkdir(parents=True, exist_ok=True)
    for k in KS:
        target = out/f'k{k:+d}.json'
        if target.exists():
            continue
        system = system_for(scen, k)
        start = time.time()
        L, R, C = system.steady_state()
        info = dict(k=k, lambda_=lam(k), D=float(system.m.D[0]), steady=steady_checks(system, L, R, C))
        states, sol = system.integrate(np.array([0.0, 50.0]), rtol=1e-8, atol=1e-10, method='RK45')
        C50 = states[-1][2]
        info['t50'] = dict(method='RK45', rtol=1e-8, atol=1e-10, nfev=int(sol.nfev))
        if args.dynamic:
            (Ld, Rd, Cd), t_end, history = system.run_to_steady_state(t_chunk=1.0 / RATES['b'],
                tol=CHECKS['dynamic_change'], t_max=1e6, method='BDF', rtol=1e-10, atol=1e-14, per_field=True)
            match = [float(np.abs(Cd[q] - C[q]).max() / C[q].max()) for q in range(len(C))]
            info['dynamic'] = dict(method='BDF', rtol=1e-10, atol=1e-14, t_chunk=1.0 / RATES['b'], t_end=t_end,
                                   last_change=history[-1][1], chunks=len(history), max_mismatch=max(match),
                                   passed=bool(max(match) <= CHECKS['dynamic_match']))
        info['seconds'] = round(time.time() - start, 1)
        np.savez_compressed(out/f'k{k:+d}.npz', C_steady=C, C_t50=C50, L_steady=L, pairs=np.array(system.m.pairs))
        write_json(target, dict(status='complete', launch=launch, finished=now(), scenario=args.scenario,
                                checks=CHECKS, **info))
        print(args.scenario, k, json.dumps(info['steady']), info.get('dynamic', ''), flush=True)


# ------------------------------------------------------------------ OT (6, 7)
def condition_name(kind, value=None):
    return kind if value is None else f'{kind}_{fmt(value)}'


def parse_condition(cond, k):
    """-> (family, parameter): inv | r0 (absolute r0) | cd (absolute dis_thr) | ref; cutoff for inv/r0 is 4 lambda."""
    if cond == 'inv':
        return 'inv', None
    if cond == 'ref':
        return 'ref', REFERENCE['dis_thr']
    kind, value = cond.split('_', 1)
    if kind == 'r0':
        return 'r0', float(value) * lam(k)
    check(kind == 'cd', f'Unknown condition {cond}')
    if value.endswith('lambda'):
        return 'cd', float(value[:-6]) * lam(k)
    return 'cd', float(value)


def cd_key(dis_thr):
    """Absolute dis_thr of a cd run, rounded so that m * lambda_k values that coincide share one run."""
    return f'cd_{dis_thr:.6f}'


def run_name(k, cond, eps, rho, nitermax):
    family, value = parse_condition(cond, k)
    if family == 'cd':
        return f'{cd_key(value)}_e{fmt(eps)}_r{fmt(rho)}_n{nitermax}'       # independent of lambda
    if family == 'ref':
        return f'ref_e{fmt(eps)}_r{fmt(rho)}_n{nitermax}'
    return f'k{k:+d}_{cond}_e{fmt(eps)}_r{fmt(rho)}_n{nitermax}'


def ot_entries(p, r0, binding, cutoff, metric=None):
    """Eligible entries (d <= cutoff, nonzero sender and receiver budgets) per pair: row, col, d, pair."""
    n = p[0].size
    xy = fields.grid_points().reshape(-1, 2)
    D = cdist(xy, xy) if metric is None else metric
    rows, cols, ds, qs = [], [], [], []
    for q, (i, j) in enumerate(pairs_of(binding)):
        x, y = np.flatnonzero(p[i].ravel() > 0), np.flatnonzero(r0[j].ravel() > 0)
        d = D[np.ix_(x, y)]
        xi, yi = np.nonzero(d <= cutoff)
        rows.append(i * n + x[xi]); cols.append(j * n + y[yi]); ds.append(d[xi, yi]); qs.append(np.full(len(xi), q))
    return (np.concatenate(rows).astype(np.int64), np.concatenate(cols).astype(np.int64), np.concatenate(ds),
            np.concatenate(qs).astype(np.int64))


def received(col, T, pair, n_pairs, n):
    return np.bincount(pair * n + col % n, weights=T, minlength=n_pairs * n).reshape(n_pairs, n)


def solve_ours(scen, family, r0, cutoff, eps, rho, nitermax, stopthr):
    from research import solve_sparse
    p, r, binding = scen['p'], scen['r0'], scen['binding']
    row, col, d, pair = ot_entries(p, r, binding, cutoff)
    c = np.log1p(d / r0) if family == 'r0' else np.where(d > H, np.log(np.maximum(d, H) / H), 0.0)
    a, b = p.reshape(len(p), -1).reshape(-1), r.reshape(len(r), -1).reshape(-1)
    try:
        T, info = solve_sparse(row, col, c, a, b, eps, rho, nitermax, stopthr)
        info = dict(check_passed=True, **{f'check_{k}': v for k, v in info.items()})
    except ValueError as error:
        return None, dict(check_passed=False, check_error=str(error), entries=len(row))
    info.update(entries=len(row), forbidden_mass=float(T[d > cutoff * (1 + 1e-12)].sum()), total_T=float(T.sum()),
                frac_ligand=float(T.sum() / a.sum()), frac_receptor=float(T.sum() / b.sum()),
                niter=info.get('check_niter'), cap_reached=info.get('check_cap_reached'))
    return received(col, T, pair, len(pairs_of(binding)), p[0].size), info


def solve_commot(scen, dis_thr, eps, rho, nitermax, stopthr, square=False):
    from scipy import sparse
    from commot_baseline import official_cot, solve_counted
    from research import sparse_kkt_check, sparse_mass_check
    p, r, binding = scen['p'], scen['r0'], scen['binding']
    nl, nr = binding.shape
    n = p[0].size
    S, R = p.reshape(nl, n).T, r.reshape(nr, n).T
    A = np.where(binding > 0, 1.0, np.inf)
    xy = fields.grid_points().reshape(-1, 2)
    M = cdist(xy, xy) ** (2 if square else 1)
    cut = dis_thr ** (2 if square else 1)
    P, iterations = solve_counted(official_cot(), 'cot_sparse', S, R, A, M, cut * np.ones_like(A), eps_p=eps,
                                  rho=rho, nitermax=nitermax, stopthr=stopthr)
    check(len(iterations) > 0, 'No iteration count recorded')
    pairs = pairs_of(binding)
    rows, cols, Ts, qs = [], [], [], []
    for q, (i, j) in enumerate(pairs):
        plan = sparse.coo_matrix(P[(i, j)])
        plan.sum_duplicates()
        rows.append(i * n + plan.row); cols.append(j * n + plan.col); Ts.append(plan.data); qs.append(np.full(plan.nnz, q))
    row, col = np.concatenate(rows).astype(np.int64), np.concatenate(cols).astype(np.int64)
    T, pair = np.maximum(np.concatenate(Ts), 0), np.concatenate(qs).astype(np.int64)
    Mval = M[row % n, col % n]
    a, b = S.T.reshape(-1), R.T.reshape(-1)
    max_amount = max(a.sum(), b.sum())
    cost_scale = float(M[M <= cut].max())
    info = dict(entries=len(row), solves=len(iterations), niter=max(iterations), nitermax=nitermax,
                cap_reached=bool(max(iterations) > nitermax), cost_scale=cost_scale,
                forbidden_mass=float(T[Mval > cut * (1 + 1e-12)].sum()), total_T=float(T.sum()),
                frac_ligand=float(T.sum() / a.sum()), frac_receptor=float(T.sum() / b.sum()))
    try:
        chk = (row, col, Mval / cost_scale, T / max_amount, a / max_amount, b / max_amount, eps, rho)
        found = dict(sparse_mass_check(*chk), **sparse_kkt_check(*chk))
        info.update(check_passed=True, **{f'check_{k}': v for k, v in found.items()})
    except ValueError as error:
        info.update(check_passed=False, check_error=str(error))
    return received(col, T, pair, len(pairs), n), info


def ot(args):
    launch = launch_info()
    scen = load_scenario(args.scenarios, args.scenario)
    name = run_name(args.k, args.condition, args.eps, args.rho, args.nitermax)
    out = Path(args.output)/args.scenario
    target = out/f'{name}.json'
    check(not target.exists(), f'Output exists: {target}')
    out.mkdir(parents=True, exist_ok=True)
    family, value = parse_condition(args.condition, args.k)
    start = time.time()
    if family in ('inv', 'r0'):
        cutoff = CUTOFF_OVER_LAMBDA * lam(args.k)
        signal, info = solve_ours(scen, family, value, cutoff, args.eps, args.rho, args.nitermax, args.stopthr)
        setting = dict(cost='log1p(d/r0)' if family == 'r0' else 'max(log(d/h), 0)', r0=value, cutoff=cutoff, h=H)
    elif family == 'cd':
        check(value >= H, f'dis_thr {value} below the grid spacing: only d = 0 eligible (plan 7.2: skipped)')
        signal, info = solve_commot(scen, value, args.eps, args.rho, args.nitermax, args.stopthr)
        setting = dict(cost='d / cost_scale (COMMOT cot_sparse)', dis_thr=value)
    else:
        check((args.eps, args.rho, args.nitermax, args.stopthr) ==
              (REFERENCE['eps'], REFERENCE['rho'], REFERENCE['nitermax'], REFERENCE['stopthr']), 'ref uses plan 7.3 settings')
        signal, info = solve_commot(scen, value, args.eps, args.rho, args.nitermax, args.stopthr, square=True)
        setting = dict(cost='d^2 / cost_scale (euc_square; archive benchmark)', dis_thr=value)
    seconds = round(time.time() - start, 1)
    if signal is not None:
        np.savez_compressed(out/f'{name}.npz', signal=signal, pairs=np.array(pairs_of(scen['binding'])))
    write_json(target, dict(status='complete', launch=launch, finished=now(), scenario=args.scenario,
        family=family, condition=args.condition, k=args.k if family in ('inv', 'r0') else None,
        eps=args.eps, rho=args.rho, nitermax=args.nitermax, stopthr=args.stopthr, seconds=seconds,
        signal_saved=signal is not None, **setting, **info))
    print(name, json.dumps({k: info.get(k) for k in ['entries', 'niter', 'cap_reached', 'check_passed']}), seconds, flush=True)


# ------------------------------------------------------------------ job list (11)
def job_conditions():
    """(k, condition) for every planned OT condition; cd at m * lambda written as cd_<m>lambda."""
    out = []
    for k in KS:
        out.append((k, 'inv'))
        out += [(k, condition_name('r0', g)) for g in R0_OVER_LAMBDA]
        out += [(k, f'cd_{fmt(m)}lambda') for m in DIS_OVER_LAMBDA if m * lam(k) >= H]
        out.append((k, condition_name('cd', DIS_FIXED)))
    return out


def jobs(args):
    """Unique runs (cd shared across lambda by absolute dis_thr), in the priority order of plan section 11."""
    folder, results = Path(args.scenarios), Path(args.results)
    manifest = json.loads((folder/'manifest.json').read_text())
    gen = [s for s in manifest['scenarios'] if s.startswith('gen_')]
    arc = [s for s in manifest['scenarios'] if s.startswith('arc_')]
    python = args.python
    lines, seen = [], set()

    def add(priority, size, scenario, k, cond, eps, rho, nitermax=NITERMAX, stopthr=STOPTHR):
        name = run_name(k, cond, eps, rho, nitermax)
        if (scenario, name) in seen:
            return
        seen.add((scenario, name))
        cmd = (f'{python} -m simulation.benchmark ot --scenarios {folder} --scenario {scenario} --k {k} '
               f'--condition {cond} --eps {fmt(eps)} --rho {fmt(rho)} --nitermax {nitermax} --stopthr {stopthr:g} '
               f'--output {results/"ot"}')
        lines.append((priority, -size, f'{scenario}/{name}', str(results/'ot'/scenario/f'{name}.json'), cmd))

    sizes = {}

    def size_of(scenario, k, cond):
        family, value = parse_condition(cond, k)
        cutoff = CUTOFF_OVER_LAMBDA * lam(k) if family in ('inv', 'r0') else value
        key = (scenario, round(cutoff, 6))
        if key not in sizes:
            scen = load_scenario(folder, scenario)
            sizes[key] = len(ot_entries(scen['p'], scen['r0'], scen['binding'], cutoff)[0])
        return sizes[key]

    for scenario, k, cond in DIAGNOSTICS:                                    # 1: D1-D3 at 200000
        add(0, size_of(scenario, k, cond), scenario, k, cond, *SETTINGS[0], nitermax=200000)
    for scenario in arc:                                                     # 1: pipeline check (R)
        add(0, 0, scenario, 0, 'ref', REFERENCE['eps'], REFERENCE['rho'], REFERENCE['nitermax'], REFERENCE['stopthr'])
    for scenario in gen:                                                     # 2: eps/rho 0.2/1
        for k, cond in job_conditions():
            add(1, size_of(scenario, k, cond), scenario, k, cond, *SETTINGS[1])
    for scenario in arc:                                                     # 3: archive scenarios, lambda_A
        for k, cond in job_conditions():
            if k == 0:
                add(2, size_of(scenario, k, cond), scenario, k, cond, *SETTINGS[0])
    for scenario in gen:                                                     # 4: main, largest first
        for k, cond in job_conditions():
            add(3, size_of(scenario, k, cond), scenario, k, cond, *SETTINGS[0])
    truth_lines = [(0, 0, f'truth/{s}', str(results/'truth'/s/'k+2.json'),
                    f'{python} -m simulation.benchmark truth --scenarios {folder} --scenario {s} --output {results/"truth"}'
                    + (' --dynamic' if s.endswith('_s1') else '')) for s in gen + arc]
    lines = truth_lines + sorted(lines, key=lambda x: (x[0], x[1]))
    table = pd.DataFrame([dict(name=n, done=d, command=c) for _, _, n, d, c in lines])
    table.to_csv(args.output, sep='\t', index=False)
    print(f'{len(table)} jobs -> {args.output}')


# ------------------------------------------------------------------ checks (10)
def h5ad_matrix(raw):
    import h5py
    with h5py.File(io.BytesIO(raw), 'r') as h:
        names = [v.decode() if isinstance(v, bytes) else str(v) for v in h['var/_index'][:]]
        X = h['X']
        if isinstance(X, h5py.Dataset):
            values = X[:]
        else:
            from scipy import sparse
            shape = tuple(X.attrs['shape'])
            values = sparse.csr_matrix((X['data'][:], X['indices'][:], X['indptr'][:]), shape=shape).toarray()
    return pd.DataFrame(values, columns=names)


def compare(a, b):
    l1 = float(np.abs(a - b).sum() / np.abs(b).sum()) if np.abs(b).sum() else np.nan
    rho = float(spearmanr(a, b).statistic) if np.ptp(a) > 0 and np.ptp(b) > 0 else np.nan
    return l1, rho


def pipeline(args):
    launch = launch_info()
    zf = open_archive(args.archive)
    rows = []
    for case in STRUCTURES:
        scenario = f'arc_c{case:02d}_sim1'
        name = run_name(0, 'ref', REFERENCE['eps'], REFERENCE['rho'], REFERENCE['nitermax'])
        with np.load(Path(args.results)/'ot'/scenario/f'{name}.npz') as z:
            signal, pairs = z['signal'], z['pairs']
        theirs = h5ad_matrix(zf.read(f'pde_model_comparison/data/case_{case}/sim_1/cot/adata.h5ad'))
        for q, (i, j) in enumerate(pairs):
            l1, rho = compare(signal[q], theirs[f'receiver-L{i+1}-R{j+1}'].to_numpy(float))
            rows.append(dict(case=case, pair=f'L{i+1}-R{j+1}', rel_l1=l1, spearman=rho,
                             meets=bool(l1 <= STABILITY['rel_l1'] and rho >= STABILITY['spearman'])))
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    table = pd.DataFrame(rows)
    table.to_csv(out/'pipeline_check.csv', index=False)
    write_json(out/'pipeline_check.json', dict(status='complete', launch=launch, finished=now(),
        thresholds=dict(rel_l1=STABILITY['rel_l1'], spearman=STABILITY['spearman']), pairs=len(table),
        failing=int((~table.meets).sum()), max_rel_l1=table.rel_l1.max(), min_spearman=table.spearman.min()))
    print(table.groupby('case')[['rel_l1', 'spearman']].agg(['max', 'min']).to_string())


def stability(args):
    launch = launch_info()
    rows, totals = [], []
    for scenario, k, cond in DIAGNOSTICS:
        loaded = []
        for nitermax in (NITERMAX, 200000):
            base = Path(args.results)/'ot'/scenario/run_name(k, cond, *SETTINGS[0], nitermax)
            meta = json.loads(base.with_suffix('.json').read_text())
            with np.load(base.with_suffix('.npz')) as z:
                loaded.append((meta, z['signal']))
        (m50, s50), (m200, s200) = loaded
        t50, t200 = s50.sum(), s200.sum()
        per = [compare(s50[q], s200[q]) for q in range(len(s50))]
        for q, (l1, rho) in enumerate(per):
            rows.append(dict(scenario=scenario, k=k, condition=cond, pair=q, rel_l1=l1, spearman=rho,
                             meets=bool(l1 <= STABILITY['rel_l1'] and rho >= STABILITY['spearman'])))
        rel_total = float(abs(t50 - t200) / t200)
        totals.append(dict(scenario=scenario, k=k, condition=cond, niter_50k=m50['niter'], niter_200k=m200['niter'],
                           cap_reached_200k=m200['cap_reached'], rel_total_difference=rel_total,
                           max_rel_l1=max(l for l, _ in per), min_spearman=min(r for _, r in per),
                           stable=bool(all(r['meets'] for r in rows if r['condition'] == cond and r['k'] == k)
                                       and rel_total <= STABILITY['total'])))
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out/'stability_pairs.csv', index=False)
    verdict = pd.DataFrame(totals)
    verdict.to_csv(out/'stability.csv', index=False)
    write_json(out/'stability.json', dict(status='complete', launch=launch, finished=now(), thresholds=STABILITY,
                                          all_stable=bool(verdict.stable.all())))
    print(verdict.to_string(index=False))


# ------------------------------------------------------------------ metrics (8)
def pair_metrics(truth, signal, r0):
    """Plan section 8 / archive evaluate_result_*.py: grid points with R0_j > 0; Spearman NaN -> 0; RMSE."""
    keep = r0.ravel() > 0
    t, s = truth.ravel()[keep], signal.ravel()[keep]
    rho = spearmanr(t, s).statistic
    nan = bool(np.isnan(rho))
    return dict(spearman=0.0 if nan else float(rho), spearman_nan=nan, rmse=float(np.sqrt(np.mean((t - s) ** 2))))


def evaluate(args):
    launch = launch_info()
    results = Path(args.results)
    rows = []
    for meta_path in sorted((results/'ot').glob('*/*.json')):
        meta = json.loads(meta_path.read_text())
        scen = load_scenario(args.scenarios, meta['scenario'])
        pairs = pairs_of(scen['binding'])
        signal = None
        if meta.get('signal_saved'):
            with np.load(meta_path.with_suffix('.npz')) as z:
                signal = z['signal']
        ks = [meta['k']] if meta['family'] in ('inv', 'r0') else list(KS)   # cd / ref shared across lambda
        for k in ks:
            tpath = results/'truth'/meta['scenario']/f'k{k:+d}.npz'
            if not tpath.exists():
                continue
            with np.load(tpath) as z:
                truths = dict(steady=z['C_steady'], t50=z['C_t50'])
            for which, C in truths.items():
                for q, (i, j) in enumerate(pairs):
                    base = dict(scenario=meta['scenario'], structure=int(scen['structure']),
                                repeat=int(scen['repeat']), source=str(scen['source']), k=k,
                                family=meta['family'], condition=meta['condition'], r0=meta.get('r0'),
                                dis_thr=meta.get('dis_thr'), eps=meta['eps'], rho=meta['rho'],
                                nitermax=meta['nitermax'], truth=which, pair=q,
                                check_passed=meta.get('check_passed'), kkt_tol=meta.get('check_kkt_tol'),
                                forbidden_mass=meta.get('forbidden_mass'), cap_reached=meta.get('cap_reached'))
                    if signal is None:
                        rows.append(dict(base, spearman=np.nan, spearman_nan=None, rmse=np.nan))
                    else:
                        rows.append(dict(base, **pair_metrics(C[q], signal[q], scen['r0'][j])))
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out/'metrics.csv.gz', index=False)
    write_json(out/'metrics.json', dict(status='complete', launch=launch, finished=now(), rows=len(rows)))


def analyze(args):
    """Plan section 9 (simulation/analysis.py) on the metrics of `evaluate`."""
    from simulation import analysis
    launch = launch_info()
    summaries = analysis.run(Path(args.metrics)/'metrics.csv.gz', args.output)
    write_json(Path(args.output)/'manifest.json', dict(status='complete', launch=launch, finished=now(),
                                                       metrics=str(args.metrics)))
    print(json.dumps(summaries, indent=2, default=str))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    c = sub.add_parser('calibrate'); c.add_argument('--archive', required=True); c.add_argument('--output', required=True)
    c = sub.add_parser('scenarios'); c.add_argument('--archive', required=True)
    c.add_argument('--calibration', required=True); c.add_argument('--output', required=True)
    c = sub.add_parser('truth'); c.add_argument('--scenarios', required=True); c.add_argument('--scenario', required=True)
    c.add_argument('--output', required=True); c.add_argument('--dynamic', action='store_true')
    c = sub.add_parser('ot'); c.add_argument('--scenarios', required=True); c.add_argument('--scenario', required=True)
    c.add_argument('--k', type=int, required=True); c.add_argument('--condition', required=True)
    c.add_argument('--eps', type=float, required=True); c.add_argument('--rho', type=float, required=True)
    c.add_argument('--nitermax', type=int, default=NITERMAX); c.add_argument('--stopthr', type=float, default=STOPTHR)
    c.add_argument('--output', required=True)
    c = sub.add_parser('jobs'); c.add_argument('--scenarios', required=True); c.add_argument('--results', required=True)
    c.add_argument('--output', required=True)
    c.add_argument('--python', default='/home/liu/.conda/envs/commot-research/bin/python')
    c = sub.add_parser('pipeline'); c.add_argument('--archive', required=True); c.add_argument('--results', required=True)
    c.add_argument('--output', required=True)
    c = sub.add_parser('stability'); c.add_argument('--results', required=True); c.add_argument('--output', required=True)
    c = sub.add_parser('evaluate'); c.add_argument('--scenarios', required=True); c.add_argument('--results', required=True)
    c.add_argument('--output', required=True)
    c = sub.add_parser('analyze'); c.add_argument('--metrics', required=True); c.add_argument('--output', required=True)
    return p


def main():
    args = parser().parse_args()
    dict(calibrate=calibrate, scenarios=scenarios, truth=truth, ot=ot, jobs=jobs, pipeline=pipeline,
         stability=stability, evaluate=evaluate, analyze=analyze)[args.command](args)


if __name__ == '__main__':
    main()
