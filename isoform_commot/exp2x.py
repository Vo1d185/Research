"""Experiment 2, expanded (reports/stage2_ppi_03/plan.md): Exp2 on the COMMOT-paper datasets.

A line (decisions; plan section 3) - Exp2's protocol unchanged except beta = 1 for every S:
  cases        case list of one dataset -> cases.json (core, length control, permutations)
  run          solve cases by name and evaluate at once (stage2_s.py functions: receptor-level received signal,
               partial Spearman against the scSeqComm target score, top 10 / 5 / 20 TFs); --save-case writes a
               stage-2 result folder for the equivalence gate
  eval-commot  receptor metrics of a commot_baseline.py output folder (reference condition, plan 3.4)
  regression   compare run outputs with Exp2's committed receptor table (post: no_ppi, raw; plan 3.6)
  summarise    per-dataset D1/D2/D3, pooled decisions over the two new datasets, references, descriptive rows

B line (descriptive; plan section 4) - the paper's protocol and the paper's metric:
  paper-ref    step 0a: the paper's metric on the paper's own received signal (paper_eval.npz of prepare-b)
  run-b        our cost on the paper panel (eps 0.1, rho 10, paper cutoff), paper metric per LR pair
  eval-b       paper metric of a commot_baseline.py output (COMMOT rerun, step 0b) and the receptor-only row;
               per-LR agreement of the rerun with the paper's received signal
  summarise-b  tables per dataset

Cost: C = log(1 + d / r0) - beta * log S_eff, r0 = 50 um, beta = 1 (beta = 0 for no_ppi).
Environments: commot-research (CPU); run-b on GPU in the torch environment (inject-gpu) with --device cuda.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import time
import zipfile

import numpy as np
import pandas as pd

from research import check
import stage2_s as ss

R0_UM = 50.0
A_CUTOFF_UM = 363.0757
STOPTHR = 1e-11
A_SETTING = dict(eps=0.2, rho=1.0, nitermax=50000)
B_SETTING = dict(eps=0.1, rho=10.0, nitermax=50000)
N_PERM = 100
TOP_TF = ss.TOP_TF                 # (10, 5, 20); first = primary
N_BOOT, BOOT_SEED, ALPHA = 1000, 1, 0.05
CORE = ['no_ppi', 'raw', 'snorm', 'uniform']
DECIDING = ('breast', 'seqfish')   # plan 3.5: decisions only on the two new datasets
REFERENCE = ('post', 'ant')        # Exp2 slices, rerun with this case list (reference, not decisions)
PAPER_DB = dict(human=('TF_TG_TRRUSTv2', 'TF_PPR_REACTOME_human'),
                mouse=('TF_TG_TRRUSTv2_RegNetwork_High_mouse', 'TF_PPR_REACTOME_mouse'))
PAPER_PPR_MIN = 1e-5
REPRO_TOL = 0.001                  # plan 4.2: |median - paper|, |mean - notebook mean|


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as s:
        for block in iter(lambda: s.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False,
                                     default=lambda x: x.item() if hasattr(x, 'item') else str(x)))


# ------------------------------------------------------------------ inputs and edges
def edges_for(spatial, cfg, cutoff):
    """Diffusion edges d <= cutoff (self included). Contact edges by config: 'hex' = stage2.spatial_edges
    (self + hexagonal neighbours, identical to Exp2); 'distance' = self + d <= contact_cutoff_um (seqFISH+);
    'none' = self only (B line: no contact pairs)."""
    import stage2
    rule = cfg.get('contact_rule', 'hex')
    if rule == 'hex':
        return dict(zip(['diffusion', 'contact'], stage2.spatial_edges(spatial, cutoff)))
    from exp2x_data import distance_pairs
    xy = spatial[['x', 'y']].to_numpy(float)
    n = len(xy)
    s, t = distance_pairs(xy, cutoff)
    s, t = np.concatenate([np.arange(n), s]), np.concatenate([np.arange(n), t])
    diffusion = (s, t, np.linalg.norm(xy[s] - xy[t], axis=1))
    if rule == 'distance':
        cs, ct = distance_pairs(xy, float(cfg['contact_cutoff_um']))
        cs, ct = np.concatenate([np.arange(n), cs]), np.concatenate([np.arange(n), ct])
    else:
        check(rule == 'none', f'Unknown contact rule {rule}')
        cs = ct = np.arange(n)
    return dict(diffusion=diffusion, contact=(cs, ct, np.linalg.norm(xy[cs] - xy[ct], axis=1)))


def load_input(source, cutoff=None):
    """stage2.load_input + edges; A line: cutoff fixed at 363.0757 um whatever config.json stores (as stage2_s);
    B line: the config's paper cutoff."""
    import stage2
    prep = stage2.verified(source)
    cfg, x, spatial, lr, ligands, receptors, a, b = stage2.load_input(source)
    stored = float(cfg['diffusion_cutoff_um'])
    cut = stored if cutoff is None else float(cutoff)
    edges = edges_for(spatial, cfg, cut)
    return prep, dict(cfg, diffusion_cutoff_um=cut, config_cutoff_um=stored), x, spatial, lr.reset_index(drop=True), \
        ligands, receptors, a, b, edges


def receptor_lengths(source, protein_map):
    """Receptor entity -> length (aa) of its longest subunit (UniProt canonical, protein_map.csv); NaN if a
    subunit is unmapped."""
    ent = pd.read_csv(Path(source)/'entities.csv', dtype=str)
    pm = pd.read_csv(protein_map, dtype={'lr_symbol': str})
    length = dict(zip(pm.lr_symbol, pm.length.astype(float)))
    return {r.entity: max((length.get(g, np.nan) for g in r.subunits.split(';')), default=np.nan)
            if all(g in length for g in r.subunits.split(';')) else np.nan for r in ent.itertuples(index=False)}


def effective(pairs, spec, rec_len=None):
    """stage2_s.effective (none | raw | snorm | uniform, permutations), plus 'length' (plan 3.3): within each
    ligand the ligand's scored S values (unscored pairs keep the type median) are reassigned so that the longest
    receptor gets the largest value; ties broken by receptor name. Distribution per ligand unchanged."""
    if spec['s'] != 'length':
        return ss.effective(pairs, spec)
    s, panel = ss.effective(pairs, dict(spec, s='raw', perm=0))
    check(panel.all(), 'length control uses the full panel')
    check(rec_len is not None, 'length control needs receptor lengths')
    scored = pairs.s_raw.notna().to_numpy()
    out = s.copy()
    for _, idx in pairs.groupby('ligand').indices.items():
        idx = [int(i) for i in idx if scored[i]]
        if len(idx) < 2:
            continue
        values = np.sort(s[idx])[::-1]
        lens = [rec_len[pairs.receptor.iloc[i]] for i in idx]
        check(np.isfinite(lens).all(), 'A scored pair has a receptor without length')
        order = sorted(range(len(idx)), key=lambda j: (-lens[j], pairs.receptor.iloc[idx[j]]))
        for rank, j in enumerate(order):
            out[idx[j]] = values[rank]
    return out, panel


# ------------------------------------------------------------------ A line: cases and run
def case_list(smoke=False):
    c = lambda name, s, beta, perm=0, group='core', keep=True: dict(
        name=name, s=s, beta=beta, perm=perm, override=None, panel='all', nitermax=A_SETTING['nitermax'],
        keep=keep, group=group)
    out = [c('no_ppi', 'none', 0.0), c('raw', 'raw', 1.0), c('snorm', 'snorm', 1.0), c('uniform', 'uniform', 1.0),
           c('length_raw', 'length', 1.0, group='descriptive')]
    for k in range(1, (2 if smoke else N_PERM) + 1):
        out += [c(f'raw_perm_{k:03d}', 'raw', 1.0, perm=k, group='perm', keep=False),
                c(f'snorm_perm_{k:03d}', 'snorm', 1.0, perm=k, group='perm', keep=False)]
    return out


def cases_cmd(args):
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    out.parent.mkdir(parents=True, exist_ok=True)
    cs = case_list(args.smoke)
    write_json(out, dict(created_utc=now(), dataset=args.dataset, smoke=args.smoke, setting=A_SETTING, r0_um=R0_UM,
                         cutoff_um=A_CUTOFF_UM, stopthr=STOPTHR, cases=cs))
    print(f'{len(cs)} cases -> {out}')


def scores_for(args, lr):
    meta = json.loads(Path(str(args.designs) + '.json').read_text())
    check(meta['input_sha256']['lr.csv'] == digest(Path(args.input)/'lr.csv'), 'Designs were made for another input')
    return ss.pair_table(lr, pd.read_csv(args.designs), pd.read_csv(args.ppi_scores), pd.read_csv(args.protein_map))


def run(args):
    import plans
    import stage2
    cases = {c['name']: c for c in json.loads(Path(args.cases).read_text())['cases']}
    todo = [cases[n] for n in args.case]
    prep, cfg, x, spatial, lr, ligands, receptors, a, b, edges = load_input(args.input, A_CUTOFF_UM)
    pairs = scores_for(args, lr)
    rec_len = receptor_lengths(args.input, args.protein_map)
    ev = ss.Evaluation(args.input, args.genes, args.db)
    n = len(x)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    st = dict(eps=A_SETTING['eps'], rho=A_SETTING['rho'])
    for case in todo:
        target = out/case['name']
        if (target/'case.json').exists():
            continue
        target.mkdir(exist_ok=True)
        s, panel = effective(pairs, case, rec_len)
        lr_case = lr[panel].reset_index(drop=True)
        row, col, c0, d, p = stage2.build_sparse_cost(lr_case, ligands, receptors, a, b, n, edges, R0_UM)
        c = c0 + case['beta'] * (-np.log(s))[p]
        start, error = time.time(), ''
        try:
            T, info = ss.solve(row, col, c, a, b, dict(st, eps=st['eps'], rho=st['rho']), case['nitermax'], args.device)
        except ValueError as e:
            T, info, error = None, {}, str(e)
        seconds = time.time() - start
        rec = dict(case=case, dataset=args.dataset, device=args.device, input=str(Path(args.input).resolve()),
                   designs_sha256=digest(args.designs), cutoff_um=cfg['diffusion_cutoff_um'],
                   config_cutoff_um=cfg['config_cutoff_um'], contact_rule=cfg.get('contact_rule', 'hex'),
                   entries=int(len(c)), seconds=round(seconds, 1), error=error, niter=info.get('niter'),
                   cap_reached=info.get('cap_reached'),
                   checks={k: v for k, v in info.items() if k.startswith(('mass', 'kkt', 'max', 'row', 'shape'))},
                   total_transport=float(T.sum()) if T is not None else None,
                   s_eff=dict(min=float(s.min()), max=float(s.max()), n_pairs=int(len(s))))
        if T is not None:
            plan = dict(sparse=True, ligands=np.asarray(ligands, str), receptors=np.asarray(receptors, str),
                        row=row, col=col, T=T)
            _, per_receptor = plans.received_signal(plan, lr_case, n)
            ev.metrics(per_receptor).to_csv(target/'metrics.csv', index=False)
            if case['keep']:
                names = sorted(per_receptor)
                np.savez_compressed(target/'signal.npz', receptors=np.asarray(names, str),
                                    signal=np.vstack([per_receptor[r] for r in names]))
            if args.save_case:
                ss.save_case(Path(args.save_case)/case['name'], prep, args.input, case, lr_case, row, col, c, T, d, p,
                             a, b, x, ligands, receptors, info)
        write_json(target/'case.json', rec)
        print(f'{now()} {args.dataset}/{case["name"]}: {len(c):,} entries, {seconds:.0f}s, {info.get("niter")} it'
              + (f' ERROR {error}' if error else ''), flush=True)


def eval_commot(args):
    """Receptor metrics (same Evaluation) of a commot_baseline.py output (case_000.npz) on the same input."""
    import plans
    import stage2
    folder, out = Path(args.results), Path(args.output)
    check(not out.exists(), f'{out} exists')
    manifest = json.loads((folder/'manifest.json').read_text())
    check(manifest.get('status') == 'complete', 'Incomplete COMMOT output')
    for f in ['expression.csv', 'spatial.csv', 'lr.csv']:
        check(manifest['input_sha256'][f] == stage2.digest(Path(args.input)/f), f'{f}: COMMOT run used another input')
    plan = plans.load(folder/'case_000.npz')
    used = pd.read_csv(folder/'used_lr.csv', dtype=str)
    n = len(plan['spots'])
    _, per_receptor = plans.received_signal(plan, used, n)
    ev = ss.Evaluation(args.input, args.genes, args.db)
    out.mkdir(parents=True)
    ev.metrics(per_receptor).to_csv(out/'metrics.csv', index=False)
    summary = pd.read_csv(folder/'summary.csv').iloc[0].to_dict()
    write_json(out/'case.json', dict(case=dict(name=args.name, group='reference'), error='', niter=summary.get('niter'),
               cap_reached=summary.get('cap_reached'), total_transport=summary.get('total_T'), seconds=summary.get('solve_seconds'),
               results=str(folder.resolve()), commot=dict(eps=summary.get('eps'), rho=summary.get('rho'),
               weights=summary.get('weights'), dis_thr_um=summary.get('dis_thr_um'), check_passed=summary.get('check_passed'))))
    print(f'{args.name}: {len(per_receptor)} receptors -> {out}')


def regression(args):
    """Plan 3.6: this script on the Exp2 posterior input must reproduce Exp2's committed receptor values."""
    ref = pd.read_csv(args.reference, index_col=0)
    res = []
    for name in args.case:
        m = pd.read_csv(Path(args.runs)/name/'metrics.csv')
        ours = m[m.top_tf.astype(str) == str(TOP_TF[0])].set_index('receptor').partial_spearman
        theirs = ref[name].dropna()
        both = ours.dropna()
        check(set(both.index) == set(theirs.index), f'{name}: receptor sets differ ({len(both)} vs {len(theirs)})')
        diff = float((both - theirs.loc[both.index]).abs().max())
        res.append(dict(case=name, receptors=int(len(both)), max_abs_diff=diff))
    passed = all(r['max_abs_diff'] <= 1e-12 for r in res)
    write_json(args.output, dict(created_utc=now(), reference=str(Path(args.reference).resolve()), cases=res,
                                 tolerance=1e-12, passed=passed))
    print(json.dumps(res), 'passed' if passed else 'FAILED')
    check(passed, 'Regression against Exp2 failed')


# ------------------------------------------------------------------ A line: decisions
def load_dataset(folder, top_tf=TOP_TF[0]):
    """Receptor x case table of partial Spearman (cases.json cases plus reference folders) and case info."""
    wide, info = ss.load_runs(folder, top_tf)
    for ref in sorted(Path(folder).glob('ref_*')):
        rec = json.loads((ref/'case.json').read_text())
        m = pd.read_csv(ref/'metrics.csv')
        wide[ref.name] = m[m.top_tf.astype(str) == str(top_tf)].set_index('receptor').partial_spearman
        info[ref.name] = dict(rec['case'], error='', niter=rec.get('niter'), cap_reached=rec.get('cap_reached'),
                              total_transport=rec.get('total_transport'), seconds=rec.get('seconds'))
    return wide, info


def common_table(wide, info):
    """Receptors defined in every core, permutation and length case (plan 3.5); references not required."""
    names = [n for n, c in info.items() if c['group'] in ('core', 'perm', 'descriptive') and n in wide]
    return wide[names].dropna()


def dataset_decisions(wide, info):
    common = common_table(wide, info)
    core_failed = [n for n in CORE + ['length_raw'] if n in info and info[n]['error']]
    M = common.mean()
    out = dict(receptors_common=int(len(common)), M={k: float(M[k]) for k in CORE + ['length_raw'] if k in M},
               failed_cases=[n for n, c in info.items() if c['error']], undetermined=bool(
                   [n for n in CORE if n in core_failed]))
    for s in ('raw', 'snorm'):
        perms = sorted(k for k in M.index if k.startswith(f'{s}_perm_'))
        if perms and s in M:
            mp = M[perms].to_numpy()
            p = float((1 + (mp >= M[s] - 1e-12).sum()) / (1 + len(mp)))
            out[f'D1_{s}'] = dict(M=float(M[s]), n_perm=len(mp), perm_mean=float(mp.mean()),
                                  perm_sd=float(mp.std(ddof=1)) if len(mp) > 1 else None, p=p,
                                  label='information' if p <= ALPHA else 'no_information')
        out[f'D2_{s}'] = ss.paired(common, s, 'no_ppi')
    out['D3_raw_vs_snorm'] = ss.paired(common, 'raw', 'snorm')
    desc = {'uniform_vs_no_ppi': ('uniform', 'no_ppi'), 'length_raw_vs_no_ppi': ('length_raw', 'no_ppi'),
            'raw_vs_length_raw': ('raw', 'length_raw')}
    out['descriptive'] = {k: ss.paired(common, a, b) for k, (a, b) in desc.items()}
    refs = [n for n, c in info.items() if c['group'] == 'reference' and n in wide]
    for r in refs:   # reference on the receptors where both are defined
        both = wide[['no_ppi', r]].dropna()
        out['descriptive'][f'no_ppi_vs_{r}'] = ss.paired(both, 'no_ppi', r)
        out['descriptive'][f'M_{r}_on_its_receptors'] = dict(M=float(both[r].mean()), M_no_ppi=float(both.no_ppi.mean()),
                                                              receptors=int(len(both)))
    return out, common


def pooled_d1(results, s):
    """Plan 3.5: Z = sum_d (M_d - mean_perm_d) / sd_perm_d, permutation seeds paired across datasets."""
    per = {}
    for ds, (wide, info, common) in results.items():
        M = common.mean()
        perms = {int(k.rsplit('_', 1)[1]): float(M[k]) for k in M.index if k.startswith(f'{s}_perm_')}
        per[ds] = (float(M[s]), perms)
    seeds = sorted(set.intersection(*[set(p) for _, p in per.values()]))
    check(len(seeds) > 1, 'Too few common permutation seeds')
    z_obs, z_perm = 0.0, np.zeros(len(seeds))
    for ds, (m, perms) in per.items():
        v = np.array([perms[k] for k in seeds])
        mu, sd = v.mean(), v.std(ddof=1)
        z_obs += (m - mu) / sd
        z_perm += (v - mu) / sd
    p = float((1 + (z_perm >= z_obs - 1e-12).sum()) / (1 + len(seeds)))
    return dict(Z=float(z_obs), n_perm=len(seeds), perm_Z_mean=float(z_perm.mean()), perm_Z_sd=float(z_perm.std(ddof=1)),
                p=p, label='information' if p <= ALPHA else 'no_information')


def pooled_paired(commons, a, b, seed=BOOT_SEED, n_boot=N_BOOT):
    """Plan 3.5: receptor-weighted mean of the per-dataset paired differences; stratified bootstrap (receptors
    resampled within each dataset; datasets in DECIDING order). improved: lower bound > 0 and every dataset's
    mean difference > 0; worse: upper bound < 0 and every dataset's < 0; otherwise uncertain."""
    diffs = {ds: (c[a] - c[b]).dropna().to_numpy() for ds, c in commons.items() if a in c and b in c}
    if len(diffs) != len(commons) or any(not len(v) for v in diffs.values()):
        return dict(delta=None, ci95=[None, None], label='not_run')
    n = {ds: len(v) for ds, v in diffs.items()}
    total = sum(n.values())
    rng = np.random.default_rng(seed)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        boots[i] = sum(diffs[ds][rng.integers(n[ds], size=n[ds])].sum() for ds in diffs) / total
    lo, hi = (float(x) for x in np.percentile(boots, [2.5, 97.5]))
    means = {ds: float(v.mean()) for ds, v in diffs.items()}
    delta = float(sum(v.sum() for v in diffs.values()) / total)
    label = ('improved' if lo > 0 and all(m > 0 for m in means.values()) else
             'worse' if hi < 0 and all(m < 0 for m in means.values()) else 'uncertain')
    return dict(delta=delta, ci95=[lo, hi], receptors=n, per_dataset_mean=means, label=label)


def conclusion(d1, d2):
    if d1 == 'information':
        return 'S carries pairing information and improves inference under this evaluation' if d2 == 'improved' else \
            'S carries pairing information but does not (yet) improve inference'
    return 'no pairing information detectable under this evaluation'


def summarise(args):
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    out.mkdir(parents=True)
    results, per_ds = {}, {}
    for ds in args.dataset:
        folder = Path(args.runs)/f'{ds}_A'
        if not (folder/'cases.json').exists():
            continue
        wide, info = load_dataset(folder)
        wide.to_csv(out/f'{ds}_receptor_partial_spearman.csv')
        dec, common = dataset_decisions(wide, info)
        for k in TOP_TF[1:]:
            w2, i2 = load_dataset(folder, k)
            c2 = common_table(w2, i2)
            dec[f'top_tf_{k}'] = {s: ss.paired(c2, s, 'no_ppi') for s in ('raw', 'snorm')}
        per_ds[ds] = dec
        results[ds] = (wide, info, common)
    pooled = {}
    groups = dict(deciding=[d for d in DECIDING if d in results], all_four=[d for d in DECIDING + REFERENCE if d in results])
    for gname, members in groups.items():
        if not members:
            continue
        sub = {d: results[d] for d in members}
        commons = {d: results[d][2] for d in members}
        undetermined = any(per_ds[d]['undetermined'] for d in members)
        g = dict(datasets=members, undetermined=undetermined)
        for s in ('raw', 'snorm'):
            g[f'D1_{s}'] = pooled_d1(sub, s) if len(members) > 1 else per_ds[members[0]].get(f'D1_{s}')
            g[f'D2_{s}'] = pooled_paired(commons, s, 'no_ppi')
            g[f'conclusion_{s}'] = 'undetermined (a core case failed)' if undetermined else \
                conclusion(g[f'D1_{s}']['label'], g[f'D2_{s}']['label'])
        d3 = pooled_paired(commons, 'raw', 'snorm')
        d3['label'] = {'improved': 'raw_better', 'worse': 'snorm_better', 'uncertain': 'indistinguishable'}.get(d3['label'], d3['label'])
        g['D3_raw_vs_snorm'] = d3
        g['descriptive'] = {k: pooled_paired(commons, a, b) for k, (a, b) in
                            {'uniform_vs_no_ppi': ('uniform', 'no_ppi'), 'length_raw_vs_no_ppi': ('length_raw', 'no_ppi'),
                             'raw_vs_length_raw': ('raw', 'length_raw')}.items()}
        pooled[gname] = g
    write_json(out/'decisions.json', dict(created_utc=now(), n_boot=N_BOOT, alpha=ALPHA, deciding=list(DECIDING),
                                          reference=list(REFERENCE), per_dataset=per_ds, pooled=pooled,
                                          note='all_four is descriptive (post/ant were seen in Exp2)'))
    print(json.dumps({g: {k: v.get('label') if isinstance(v, dict) else v for k, v in d.items()
                          if k.startswith(('D1', 'D2', 'D3', 'conclusion'))} for g, d in pooled.items()}, indent=2))


# ------------------------------------------------------------------ B line: paper metric
def paper_items(lr, species, zip_path):
    """get_downstream_genes of cell-level_downstream_gene_correlation.ipynb (tf_ppr_thresh = 1e-5), on the scSeqComm
    tables exported in the paper's Zenodo archive: receptor 'A_B' -> 'A,B' matched exactly to the REACTOME PPR
    receptor names; pairs whose receptor has no TF with a target list are skipped; targets = concatenated targets
    of the receptor's listed TFs (duplicates kept; only membership is used)."""
    tg_name, ppr_name = PAPER_DB[species]
    base = 'downstream_gene_correlation/scseqcomm/data/'
    with zipfile.ZipFile(zip_path) as z:
        ppr = pd.read_csv(io.BytesIO(z.read(f'{base}{ppr_name}.csv')), index_col=0)
        ppr = ppr.loc[ppr['tf_PPR'] >= PAPER_PPR_MIN, :]
        recs = set(ppr['receptor'].values)
        avail = list(pd.read_csv(io.BytesIO(z.read(f'{base}TF_list-{tg_name}.txt')))['x'].values)
        cache, items = {}, []
        for lig, rec in zip(lr.ligand, lr.receptor):
            name = ','.join(rec.split('_'))
            if name not in recs:
                continue
            tfs = list(ppr.loc[ppr['receptor'] == name, 'tf'].values)
            if not set(tfs) & set(avail):
                continue
            targets = []
            for tf in tfs:
                if tf in avail:
                    if tf not in cache:
                        cache[tf] = list(pd.read_csv(io.BytesIO(z.read(f'{base}{tg_name}/{tf}.csv')))['x'].values)
                    targets.extend(cache[tf])
            items.append((lig, rec, targets))
    return items


def paper_metric(signal, E, genes, items):
    """cell_tg_ranking_analysis: per LR pair, Spearman over all spots between the received signal and the mean
    expression of its targets among the evaluation genes (HVG); NaN signal/target -> 0; NaN correlation -> 0."""
    from scipy.stats import spearmanr
    genes = pd.Index(genes)
    rows = []
    n = E.shape[0]
    for lig, rec, targets in items:
        sig = signal.get((lig, rec))
        sig = np.zeros(n) if sig is None else np.nan_to_num(np.asarray(sig, float))
        mask = genes.isin(targets)
        tgt = np.nan_to_num(E[:, mask].mean(axis=1)) if mask.any() else np.zeros(n)
        with np.errstate(all='ignore'):
            c = spearmanr(sig, tgt)[0] if np.ptp(sig) > 0 and np.ptp(tgt) > 0 else np.nan
        rows.append(dict(ligand=lig, receptor=rec, n_targets_eval=int(mask.sum()), signal_zero=bool(not sig.any()),
                         spearman=0.0 if np.isnan(c) else float(c)))
    return pd.DataFrame(rows)


def summary_row(m):
    return dict(n=int(len(m)), mean=float(m.spearman.mean()), median=float(m.spearman.median()))


def load_paper_eval(inp):
    import stage2
    prep = stage2.verified(inp)
    check(digest(Path(inp)/'paper_eval.npz') == prep['paper_eval']['sha256'], 'paper_eval.npz changed')
    z = np.load(Path(inp)/'paper_eval.npz', allow_pickle=False)
    lr = pd.read_csv(Path(inp)/'lr.csv', dtype=str)
    return prep, z, lr


def paper_ref(args):
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    prep, z, lr = load_paper_eval(args.input)
    p = prep['paper']
    items = paper_items(lr, p['species'], args.zenodo_zip)
    sig = {(l, r): z['paper_signal'][i] for i, (l, r) in enumerate(zip(z['ligand'], z['receptor']))}
    m = paper_metric(sig, z['E'], z['genes'], items)
    out.mkdir(parents=True)
    m.to_csv(out/'paper_output_per_lr.csv', index=False)
    s = summary_row(m)
    checks = dict(n_eval=s['n'] == p['n_eval'], median=abs(s['median'] - p['median']) <= REPRO_TOL,
                  mean=abs(s['mean'] - p['mean']) <= REPRO_TOL, n_hvg=prep['paper_eval']['n_hvg_matches'])
    write_json(out/'paper_ref.json', dict(created_utc=now(), dataset=prep['dataset'], ours=s,
               paper=dict(n_eval=p['n_eval'], median=p['median'], mean_notebook=p['mean'], n_hvg=p['n_hvg']),
               n_hvg=prep['paper_eval']['n_hvg'], tolerance=REPRO_TOL, checks=checks, reproduced=all(checks.values())))
    print(f'{prep["dataset"]} 0a: n={s["n"]} mean={s["mean"]:.5f} median={s["median"]:.5f} '
          f'(paper n={p["n_eval"]} mean={p["mean"]:.5f} median={p["median"]:.3f}) '
          + ('reproduced' if all(checks.values()) else f'NOT reproduced {checks}'), flush=True)


def run_b(args):
    import plans
    import stage2
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    prep, cfg, x, spatial, lr, ligands, receptors, a, b, edges = load_input(args.input)
    check(cfg.get('line') == 'B' and (lr.type == 'diffusion').all(), 'run-b needs a prepare-b input')
    _, z, _ = load_paper_eval(args.input)
    items = paper_items(lr, prep['paper']['species'], args.zenodo_zip)
    pairs = scores_for(args, lr)
    n = len(x)
    names = dict(ours_no_ppi=('none', 0.0), ours_raw=('raw', 1.0), ours_snorm=('snorm', 1.0), ours_uniform=('uniform', 1.0))
    for name in args.case:
        target = out/name
        if (target/'case.json').exists():
            continue
        target.mkdir(exist_ok=True)
        s_kind, beta = names[name]
        s, panel = ss.effective(pairs, dict(s=s_kind, perm=0, override=None, panel='all'))
        row, col, c0, d, p = stage2.build_sparse_cost(lr, ligands, receptors, a, b, n, edges, R0_UM)
        c = c0 + beta * (-np.log(s))[p]
        start, error = time.time(), ''
        try:
            T, info = ss.solve(row, col, c, a, b, dict(eps=B_SETTING['eps'], rho=B_SETTING['rho']),
                               args.nitermax or B_SETTING['nitermax'], args.device)
        except ValueError as e:
            T, info, error = None, {}, str(e)
        seconds = time.time() - start
        if T is not None:
            plan = dict(sparse=True, ligands=np.asarray(ligands, str), receptors=np.asarray(receptors, str),
                        row=row, col=col, T=T)
            per_pair, _ = plans.received_signal(plan, lr, n)
            paper_metric(per_pair, z['E'], z['genes'], items).to_csv(target/'per_lr.csv', index=False)
            np.savez_compressed(target/'signal.npz', ligand=lr.ligand.to_numpy(str), receptor=lr.receptor.to_numpy(str),
                                signal=np.vstack([per_pair[(l, r)] for l, r in zip(lr.ligand, lr.receptor)]))
        write_json(target/'case.json', dict(case=name, s=s_kind, beta=beta, dataset=prep['dataset'], device=args.device,
                   eps=B_SETTING['eps'], rho=B_SETTING['rho'], nitermax=args.nitermax or B_SETTING['nitermax'],
                   cutoff_um=cfg['diffusion_cutoff_um'], r0_um=R0_UM, entries=int(len(c)), seconds=round(seconds, 1),
                   error=error, niter=info.get('niter'), cap_reached=info.get('cap_reached'),
                   total_transport=float(T.sum()) if T is not None else None, designs_sha256=digest(args.designs)))
        print(f'{now()} B {prep["dataset"]}/{name}: {len(c):,} entries, {seconds:.0f}s, {info.get("niter")} it'
              + (f' ERROR {error}' if error else ''), flush=True)


def eval_b(args):
    """Step 0b and the receptor-only row: paper metric of the COMMOT rerun (commot_baseline.py output), per-LR
    agreement with the paper's received signal (Spearman over spots; relative L1), and receptor expression as signal."""
    import plans
    import stage2
    from scipy.stats import spearmanr
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    prep, z, lr = load_paper_eval(args.input)
    items = paper_items(lr, prep['paper']['species'], args.zenodo_zip)
    folder = Path(args.results)
    manifest = json.loads((folder/'manifest.json').read_text())
    for f in ['expression.csv', 'spatial.csv', 'lr.csv']:
        check(manifest['input_sha256'][f] == stage2.digest(Path(args.input)/f), f'{f}: COMMOT run used another input')
    plan = plans.load(folder/'case_000.npz')
    used = pd.read_csv(folder/'used_lr.csv', dtype=str)
    n = len(plan['spots'])
    per_pair, _ = plans.received_signal(plan, used, n)
    out.mkdir(parents=True)
    m = paper_metric(per_pair, z['E'], z['genes'], items)
    m.to_csv(out/'commot_rerun_per_lr.csv', index=False)
    expr = pd.read_csv(Path(args.input)/'expression.csv', index_col=0)
    rec_only = paper_metric({(l, r): expr[r].to_numpy(float) for l, r in zip(lr.ligand, lr.receptor)},
                            z['E'], z['genes'], items)
    rec_only.to_csv(out/'receptor_only_per_lr.csv', index=False)
    agree = []
    for i, (l, r) in enumerate(zip(z['ligand'], z['receptor'])):
        ref = z['paper_signal'][i]
        ours = per_pair.get((l, r), np.zeros(n))
        with np.errstate(all='ignore'):
            rho = spearmanr(ours, ref)[0] if np.ptp(ours) > 0 and np.ptp(ref) > 0 else np.nan
        agree.append(dict(ligand=l, receptor=r, spearman_vs_paper=rho,
                          rel_l1=float(np.abs(ours - ref).sum() / max(np.abs(ref).sum(), 1e-300)),
                          both_zero=bool(not ours.any() and not ref.any())))
    agree = pd.DataFrame(agree)
    agree.to_csv(out/'rerun_vs_paper_signal.csv', index=False)
    summary = pd.read_csv(folder/'summary.csv').iloc[0].to_dict()
    write_json(out/'eval_b.json', dict(created_utc=now(), dataset=prep['dataset'], commot_rerun=summary_row(m),
               receptor_only=summary_row(rec_only), commot_summary={k: summary.get(k) for k in
               ['eps', 'rho', 'weights', 'dis_thr_um', 'niter', 'solves', 'solves_at_cap', 'solve_seconds']},
               agreement=dict(pairs=int(len(agree)), defined=int(agree.spearman_vs_paper.notna().sum()),
                              both_zero=int(agree.both_zero.sum()),
                              spearman_median=float(agree.spearman_vs_paper.median()),
                              spearman_min=float(agree.spearman_vs_paper.min()),
                              rel_l1_median=float(agree.rel_l1.median()), rel_l1_max=float(agree.rel_l1.max()),
                              reproduced=bool(agree.spearman_vs_paper.median() >= 0.99))))
    print(f'{prep["dataset"]} 0b: rerun {summary_row(m)}, receptor-only {summary_row(rec_only)}, '
          f'signal Spearman vs paper median {agree.spearman_vs_paper.median():.4f}', flush=True)


def paired_lr(a, b, seed=BOOT_SEED, n_boot=N_BOOT):
    d = (a.set_index(['ligand', 'receptor']).spearman - b.set_index(['ligand', 'receptor']).spearman).dropna().to_numpy()
    rng = np.random.default_rng(seed)
    boots = d[rng.integers(len(d), size=(n_boot, len(d)))].mean(1)
    lo, hi = (float(x) for x in np.percentile(boots, [2.5, 97.5]))
    return dict(delta_mean=float(d.mean()), ci95=[lo, hi], delta_median=float(np.median(d)), pairs=int(len(d)))


def summarise_b(args):
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    out.mkdir(parents=True)
    res = {}
    for ds in args.dataset:
        root = Path(args.runs)
        tables = {}
        ref = root/f'{ds}_B_paper_ref'/'paper_output_per_lr.csv'
        if ref.exists():
            tables['paper_output'] = pd.read_csv(ref)
        ev = root/f'{ds}_B_eval'
        for f, name in [('commot_rerun_per_lr.csv', 'commot_rerun'), ('receptor_only_per_lr.csv', 'receptor_only')]:
            if (ev/f).exists():
                tables[name] = pd.read_csv(ev/f)
        for name in ['ours_no_ppi', 'ours_raw', 'ours_snorm', 'ours_uniform']:
            f = root/f'{ds}_B'/name/'per_lr.csv'
            if f.exists():
                tables[name] = pd.read_csv(f)
        if not tables:
            continue
        rows = {k: summary_row(v) for k, v in tables.items()}
        comp = {}
        for a_ in ['ours_raw', 'ours_snorm', 'ours_uniform']:
            if a_ in tables and 'ours_no_ppi' in tables:
                comp[f'{a_}_vs_ours_no_ppi'] = paired_lr(tables[a_], tables['ours_no_ppi'])
        for a_ in ['ours_no_ppi', 'ours_raw', 'ours_snorm', 'ours_uniform', 'receptor_only']:
            if a_ in tables and 'commot_rerun' in tables:
                comp[f'{a_}_vs_commot_rerun'] = paired_lr(tables[a_], tables['commot_rerun'])
        if 'commot_rerun' in tables and 'paper_output' in tables:
            comp['commot_rerun_vs_paper_output'] = paired_lr(tables['commot_rerun'], tables['paper_output'])
        repro = {}
        for f, k in [(root/f'{ds}_B_paper_ref'/'paper_ref.json', '0a'), (root/f'{ds}_B_eval'/'eval_b.json', '0b')]:
            if f.exists():
                repro[k] = json.loads(f.read_text())
        res[ds] = dict(metric=rows, comparisons=comp, reproduction=dict(
            step0a=repro.get('0a', {}).get('reproduced'), step0b=repro.get('0b', {}).get('agreement', {}).get('reproduced'),
            step0a_detail=repro.get('0a'), step0b_detail=repro.get('0b', {}).get('agreement')))
        pd.concat([v.assign(condition=k) for k, v in tables.items()]).to_csv(out/f'{ds}_per_lr.csv', index=False)
    write_json(out/'b_line.json', dict(created_utc=now(), note='descriptive only (plan section 4)', datasets=res))
    print(json.dumps({d: v['metric'] for d, v in res.items()}, indent=2))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    c = sub.add_parser('cases'); c.add_argument('--dataset', required=True); c.add_argument('--output', required=True)
    c.add_argument('--smoke', action='store_true')
    r = sub.add_parser('run')
    for k in ('input', 'designs', 'ppi-scores', 'protein-map', 'genes', 'db', 'cases', 'output', 'dataset'):
        r.add_argument('--' + k, required=True)
    r.add_argument('--case', nargs='+', required=True); r.add_argument('--device', default='cpu')
    r.add_argument('--save-case')
    e = sub.add_parser('eval-commot')
    for k in ('input', 'genes', 'db', 'results', 'output', 'name'):
        e.add_argument('--' + k, required=True)
    g = sub.add_parser('regression')
    for k in ('runs', 'reference', 'output'):
        g.add_argument('--' + k, required=True)
    g.add_argument('--case', nargs='+', required=True)
    s = sub.add_parser('summarise'); s.add_argument('--runs', required=True); s.add_argument('--output', required=True)
    s.add_argument('--dataset', nargs='+', default=list(DECIDING + REFERENCE))
    q = sub.add_parser('paper-ref')
    for k in ('input', 'zenodo-zip', 'output'):
        q.add_argument('--' + k, required=True)
    b = sub.add_parser('run-b')
    for k in ('input', 'designs', 'ppi-scores', 'protein-map', 'zenodo-zip', 'output'):
        b.add_argument('--' + k, required=True)
    b.add_argument('--case', nargs='+', required=True); b.add_argument('--device', default='cuda')
    b.add_argument('--nitermax', type=int, help='smoke only; default 50000')
    v = sub.add_parser('eval-b')
    for k in ('input', 'zenodo-zip', 'results', 'output'):
        v.add_argument('--' + k, required=True)
    t = sub.add_parser('summarise-b'); t.add_argument('--runs', required=True); t.add_argument('--output', required=True)
    t.add_argument('--dataset', nargs='+', default=['post', 'breast', 'seqfish'])
    return p


if __name__ == '__main__':
    a = parser().parse_args()
    dict(cases=cases_cmd, run=run, **{'eval-commot': eval_commot}, regression=regression, summarise=summarise,
         **{'paper-ref': paper_ref, 'run-b': run_b, 'eval-b': eval_b, 'summarise-b': summarise_b})[a.command](a)
