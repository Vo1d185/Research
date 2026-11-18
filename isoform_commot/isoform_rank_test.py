"""Isoform rank test (reports/isoform_rank_test/plan.md): PLM-interact score of one ligand against two
isoforms of the same receptor gene, normalised within the group: S' = S / (S of both rows of the group).

  fetch     UniProt: ligand canonical sequence (entry JSON), receptor isoform sequence by its own isoform id
            (FASTA of e.g. P21802-3), plus the receptor entries' ALTERNATIVE PRODUCTS -> data/isoform_rank_test/uniprot/
  build     sequences.csv (versions, lengths, isoform names), receptor_check.csv (IIIb vs IIIc sequences
            must differ) and jobs.csv (each pair in both orders) -> reports/isoform_rank_test/
  evaluate  plm_score.py outputs of each model -> rows.csv, groups.csv, summary.csv, tables.md
"""
import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from math import comb
from pathlib import Path
from urllib.request import Request, urlopen

import pandas as pd

from research import check, launch_info

UA = 'isoform-commot-research/1.0'
REQUIRED = ['group_id', 'ligand_gene', 'ligand_acc', 'receptor_gene', 'receptor_acc', 'receptor_isoform',
            'isoform_uniprot_id_to_confirm', 'expected', 'confidence']
EXPECTED = {'bind', 'nobind', 'weak', 'control_both'}
MAX_TOKENS = 1603


def sha256(text):
    return hashlib.sha256(text.encode()).hexdigest()


def get(url):
    with urlopen(Request(url, headers={'User-Agent': UA}), timeout=60) as r:
        return r.read().decode(), {k.lower(): v for k, v in r.headers.items()}


def read_testset(path):
    t = pd.read_csv(path, dtype=str, keep_default_na=False)
    for c in REQUIRED:
        check(c in t, f'testset lacks column {c}')
        empty = t.index[t[c].str.strip() == ''].tolist()
        check(not empty, f'empty {c} in rows {[i + 2 for i in empty]} (file line numbers); ask the user')
    t = t.apply(lambda s: s.str.strip())
    check(set(t.expected) <= EXPECTED, f'unknown expected values {set(t.expected) - EXPECTED}')
    check(set(t.receptor_isoform) <= {'IIIb', 'IIIc'}, 'receptor_isoform must be IIIb or IIIc')
    for r in t.itertuples():
        check(re.fullmatch(re.escape(r.receptor_acc) + r'-\d+', r.isoform_uniprot_id_to_confirm),
              f'{r.isoform_uniprot_id_to_confirm} is not an isoform id of {r.receptor_acc}')
    for g, rows in t.groupby('group_id'):
        check(len(rows) == 2 and set(rows.receptor_isoform) == {'IIIb', 'IIIc'},
              f'{g}: a group must have one IIIb and one IIIc row')
        check(rows.ligand_acc.nunique() == 1 and rows.receptor_acc.nunique() == 1, f'{g}: mixed ligand or receptor')
    iso = t.groupby(['receptor_acc', 'receptor_isoform']).isoform_uniprot_id_to_confirm.nunique()
    check((iso == 1).all(), 'one receptor isoform label must map to one isoform id')
    return t


def kind(rows):
    e = set(rows.expected)
    if e == {'bind', 'nobind'}:
        return 'directional'
    if e == {'control_both'}:
        return 'control'
    check('weak' in e, f'group with expected {e} is neither directional, control nor weak')
    return 'weak'


# ---------------------------------------------------------------- fetch
def fetch(args):
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    t = read_testset(args.testset)
    out.mkdir(parents=True)
    log, releases = [], set()
    for acc in sorted(set(t.ligand_acc) | set(t.receptor_acc)):
        url = f'https://rest.uniprot.org/uniprotkb/{acc}.json'
        text, h = get(url)
        (out/f'{acc}.json').write_text(text)
        releases.add(h.get('x-uniprot-release'))
        log.append(dict(id=acc, url=url, file=f'{acc}.json', sha256=sha256(text), release=h.get('x-uniprot-release')))
    for iso in sorted(set(t.isoform_uniprot_id_to_confirm)):
        url = f'https://rest.uniprot.org/uniprotkb/{iso}.fasta'
        text, h = get(url)
        (out/f'{iso}.fasta').write_text(text)
        releases.add(h.get('x-uniprot-release'))
        log.append(dict(id=iso, url=url, file=f'{iso}.fasta', sha256=sha256(text), release=h.get('x-uniprot-release')))
    check(len(releases) == 1, f'several UniProt releases during the fetch: {releases}')
    (out/'fetch.json').write_text(json.dumps(dict(fetched_utc=datetime.now(timezone.utc).isoformat(timespec='seconds'),
        uniprot_release=releases.pop(), testset=args.testset, user_agent=UA, files=log), indent=2))


# ---------------------------------------------------------------- build
def parse_fasta(text, expected_id):
    lines = text.strip().splitlines()
    check(len(lines) >= 2 and lines[0].startswith('>'), f'{expected_id}: not a FASTA record')
    check(sum(l.startswith('>') for l in lines) == 1, f'{expected_id}: several FASTA records')
    check(f'|{expected_id}|' in lines[0], f'{expected_id}: header is {lines[0]}')
    return lines[0], ''.join(lines[1:]).strip().upper()


def alt_products(entry):
    """isoform id -> (name, synonyms, sequence status) from the ALTERNATIVE PRODUCTS comment."""
    out = {}
    for c in entry.get('comments', []):
        if c.get('commentType') != 'ALTERNATIVE PRODUCTS':
            continue
        for i in c.get('isoforms', []):
            name = i.get('name', {}).get('value', '')
            syn = ';'.join(s.get('value', '') for s in i.get('synonyms', []))
            for x in i.get('isoformIds', []):
                out[x] = (name, syn, i.get('isoformSequenceStatus', ''))
    return out


def build(args):
    src, out = Path(args.uniprot), Path(args.output)
    t = read_testset(args.testset)
    fetched = json.loads((src/'fetch.json').read_text())
    entries = {a: json.loads((src/f'{a}.json').read_text()) for a in set(t.ligand_acc) | set(t.receptor_acc)}
    rows = []
    for acc in sorted(set(t.ligand_acc)):
        e = entries[acc]
        check(e['primaryAccession'] == acc, f'{acc}: entry is {e["primaryAccession"]} (merged or demerged?)')
        gene = t.loc[t.ligand_acc == acc, 'ligand_gene'].iloc[0]
        rows.append(dict(id=acc, role='ligand', gene=gene, label='canonical', uniprot_isoform_name='canonical',
                         uniprot_synonyms='', sequence_status='canonical (entry sequence)',
                         entry_version=e['entryAudit']['entryVersion'], sequence_version=e['entryAudit']['sequenceVersion'],
                         sequence=e['sequence']['value'].upper()))
    for (acc, gene, label, iso), _ in t.groupby(['receptor_acc', 'receptor_gene', 'receptor_isoform',
                                                 'isoform_uniprot_id_to_confirm']):
        e = entries[acc]
        check(e['primaryAccession'] == acc, f'{acc}: entry is {e["primaryAccession"]}')
        header, seq = parse_fasta((src/f'{iso}.fasta').read_text(), iso)
        ap = alt_products(e)
        check(iso in ap, f'{iso} is not listed in the ALTERNATIVE PRODUCTS of {acc}')
        name, syn, status = ap[iso]
        canonical = e['sequence']['value'].upper()
        check((seq == canonical) == (status == 'Displayed'),
              f'{iso}: status {status} but sequence {"equals" if seq == canonical else "differs from"} the canonical')
        rows.append(dict(id=iso, role='receptor', gene=gene, label=label, uniprot_isoform_name=name,
                         uniprot_synonyms=syn, sequence_status=status + (' (= canonical)' if seq == canonical else ''),
                         entry_version=e['entryAudit']['entryVersion'], sequence_version=e['entryAudit']['sequenceVersion'],
                         sequence=seq))
    seqs = pd.DataFrame(rows)
    check(set(''.join(seqs.sequence)) <= set('ACDEFGHIKLMNPQRSTVWYXBUZO'), 'letters outside the ESM vocabulary')
    seqs['length'] = seqs.sequence.str.len()
    seqs['sequence_sha256'] = seqs.sequence.map(sha256)
    seqs['uniprot_release'] = fetched['uniprot_release']
    by_id = seqs.set_index('id')
    checks = []
    for gene, g in seqs[seqs.role == 'receptor'].groupby('gene'):
        b, c = g[g.label == 'IIIb'].iloc[0], g[g.label == 'IIIc'].iloc[0]
        check(b.sequence != c.sequence, f'{gene}: IIIb and IIIc sequences are identical')
        pre = next((i for i, (x, y) in enumerate(zip(b.sequence, c.sequence)) if x != y), min(b.length, c.length))
        suf = next((i for i, (x, y) in enumerate(zip(b.sequence[::-1], c.sequence[::-1])) if x != y),
                   min(b.length, c.length))
        checks.append(dict(gene=gene, IIIb_id=b.id, IIIb_length=b.length, IIIc_id=c.id, IIIc_length=c.length,
                           identical=False, length_difference=int(b.length - c.length), common_prefix=pre,
                           common_suffix=suf, differing_region_IIIb=f'{pre + 1}-{b.length - suf}',
                           differing_region_IIIc=f'{pre + 1}-{c.length - suf}'))
    jobs = []
    for r in t.drop_duplicates(['ligand_acc', 'isoform_uniprot_id_to_confirm']).itertuples():
        lig, rec = r.ligand_acc, r.isoform_uniprot_id_to_confirm
        for a, b in [(lig, rec), (rec, lig)]:
            jobs.append(dict(pair_id=f'{a}__{b}', query=by_id.loc[a, 'sequence'], text=by_id.loc[b, 'sequence']))
    jobs = pd.DataFrame(jobs)
    n_tokens = jobs['query'].str.len() + jobs['text'].str.len() + 3
    out.mkdir(parents=True, exist_ok=True)
    for name in ['sequences.csv', 'receptor_check.csv', 'jobs.csv', 'build.json']:
        check(not (out/name).exists(), f'{out/name} exists')
    seqs.to_csv(out/'sequences.csv', index=False)
    pd.DataFrame(checks).to_csv(out/'receptor_check.csv', index=False)
    jobs.to_csv(out/'jobs.csv', index=False)
    groups = {g: kind(rows) for g, rows in t.groupby('group_id')}
    (out/'build.json').write_text(json.dumps(dict(status='complete', launch=launch_info(), testset=args.testset,
        testset_sha256=hashlib.sha256(Path(args.testset).read_bytes()).hexdigest(), uniprot=fetched,
        rows=len(t), groups=len(groups), group_kinds=pd.Series(groups).value_counts().to_dict(),
        jobs=len(jobs), max_tokens=int(n_tokens.max()), over_official_limit=int((n_tokens > MAX_TOKENS).sum())),
        indent=2, default=str))
    print(pd.DataFrame(checks).to_string(index=False))
    print(seqs.drop(columns=['sequence', 'sequence_sha256']).to_string(index=False))


# ---------------------------------------------------------------- evaluate
def binom_tail(k, n):
    """P(X >= k) for X ~ Binomial(n, 0.5): descriptive only (plan 4)."""
    return sum(comb(n, i) for i in range(k, n + 1)) / 2 ** n


def evaluate(args):
    out = Path(args.output)
    check(not out.exists(), f'{out} exists')
    t = read_testset(args.testset)
    t['kind'] = t.group_id.map({g: kind(rows) for g, rows in t.groupby('group_id')})
    rows, groups, summary = [], [], []
    for spec in args.scores:
        model, path = spec.split('=', 1)
        meta = json.loads(Path(path + '.meta.json').read_text())
        check(meta['status'] == 'complete', f'{path} incomplete')
        check(meta['input_sha256'] == hashlib.sha256(Path(args.jobs).read_bytes()).hexdigest(),
              f'{path} was not scored from {args.jobs}')
        s = pd.read_csv(path).set_index('pair_id').score
        check(s.notna().all() and ((s >= 0) & (s <= 1)).all(), f'{path}: scores outside [0, 1]')
        m = t.copy()
        m['model'] = model
        m['score_lr'] = [s[f'{a}__{b}'] for a, b in zip(m.ligand_acc, m.isoform_uniprot_id_to_confirm)]
        m['score_rl'] = [s[f'{b}__{a}'] for a, b in zip(m.ligand_acc, m.isoform_uniprot_id_to_confirm)]
        m['score'] = (m.score_lr + m.score_rl) / 2
        total = m.groupby('group_id').score.transform('sum')
        m['s_norm'] = m.score / total.where(total > 0)          # plan 3.1: undefined if the group sum is 0
        rows.append(m)
        for g, r in m.groupby('group_id', sort=False):
            b, c = r[r.receptor_isoform == 'IIIb'].iloc[0], r[r.receptor_isoform == 'IIIc'].iloc[0]
            row = dict(model=model, group_id=g, kind=r.kind.iloc[0], ligand=r.ligand_gene.iloc[0],
                       receptor=r.receptor_gene.iloc[0], confidence=';'.join(sorted(set(r.confidence))),
                       score_IIIb=b.score, score_IIIc=c.score, s_norm_IIIb=b.s_norm, s_norm_IIIc=c.s_norm,
                       diff_IIIb_minus_IIIc=b.s_norm - c.s_norm)
            if row['kind'] == 'directional':
                bind, nob = r[r.expected == 'bind'].iloc[0], r[r.expected == 'nobind'].iloc[0]
                row.update(expected_higher=bind.receptor_isoform, s_norm_bind=bind.s_norm, s_norm_nobind=nob.s_norm,
                           diff_bind_minus_nobind=bind.s_norm - nob.s_norm,
                           hit=bool(pd.notna(bind.s_norm) and bind.s_norm > nob.s_norm))   # plan 3.2: strict >
            if row['kind'] == 'control':
                row['abs_diff'] = abs(b.s_norm - c.s_norm)
            groups.append(row)
    rows, groups = pd.concat(rows, ignore_index=True), pd.DataFrame(groups)
    for model, g in groups.groupby('model', sort=False):
        d = g[g.kind == 'directional']
        for subset, x in [('all', d), ('high confidence', d[d.confidence == 'high'])]:
            k, n = int(x.hit.sum()), len(x)
            iiic = int((x.expected_higher == 'IIIc').sum())
            summary.append(dict(model=model, subset=subset, groups=n, hits=k, hit_rate=k / n if n else float('nan'),
                                random_baseline=0.5, binom_p_ge=binom_tail(k, n) if n else float('nan'),
                                always_IIIc_higher_rate=iiic / n if n else float('nan'),
                                always_IIIb_higher_rate=(n - iiic) / n if n else float('nan'),
                                median_diff_bind_minus_nobind=float(x.diff_bind_minus_nobind.median()) if n else float('nan')))
    summary = pd.DataFrame(summary)
    sign = (groups.groupby(['model', 'receptor'], sort=False)
            .diff_IIIb_minus_IIIc.agg(groups='size', IIIb_higher=lambda v: int((v > 0).sum()),
                                      IIIc_higher=lambda v: int((v < 0).sum())).reset_index())
    out.mkdir(parents=True)
    rows.to_csv(out/'rows.csv', index=False)
    groups.to_csv(out/'groups.csv', index=False)
    summary.to_csv(out/'summary.csv', index=False)
    sign.to_csv(out/'receptor_sign.csv', index=False)
    (out/'manifest.json').write_text(json.dumps(dict(status='complete', launch=launch_info(), testset=args.testset,
        jobs=args.jobs, scores=args.scores, created_utc=datetime.now(timezone.utc).isoformat(timespec='seconds')),
        indent=2, default=str))
    f = lambda v: '—' if pd.isna(v) else f'{v:.4f}'
    md = []
    for model in groups.model.unique():
        g = groups[groups.model == model]
        md += [f'### {model}', '', '| 组 | 类型 | 置信度 | 配体 | 受体 | S(IIIb) | S(IIIc) | S\'(IIIb) | S\'(IIIc) | 预期较高 | S\'(bind) − S\'(nobind) | 命中 |',
               '|---|---|---|---|---|---|---|---|---|---|---|---|']
        for r in g.itertuples():
            hit = ('是' if r.hit else '否') if r.kind == 'directional' else '—'
            md.append(f'| {r.group_id} | {r.kind} | {r.confidence} | {r.ligand} | {r.receptor} | {r.score_IIIb:.4g} | '
                      f'{r.score_IIIc:.4g} | {f(r.s_norm_IIIb)} | {f(r.s_norm_IIIc)} | '
                      f'{getattr(r, "expected_higher", "") if r.kind == "directional" else "—"} | '
                      f'{f(r.diff_bind_minus_nobind) if r.kind == "directional" else "—"} | {hit} |')
        md.append('')
    (out/'tables.md').write_text('\n'.join(md))
    print(summary.to_string(index=False))
    print(sign.to_string(index=False))
    print(groups[groups.kind != 'directional'].to_string(index=False))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    f = sub.add_parser('fetch')
    f.add_argument('--testset', default='reports/isoform_rank_test/testset.csv')
    f.add_argument('--output', default='data/isoform_rank_test/uniprot')
    b = sub.add_parser('build')
    b.add_argument('--testset', default='reports/isoform_rank_test/testset.csv')
    b.add_argument('--uniprot', default='data/isoform_rank_test/uniprot')
    b.add_argument('--output', default='reports/isoform_rank_test')
    e = sub.add_parser('evaluate')
    e.add_argument('--testset', default='reports/isoform_rank_test/testset.csv')
    e.add_argument('--jobs', default='reports/isoform_rank_test/jobs.csv')
    e.add_argument('--scores', nargs='+', required=True, help='model=path to a plm_score.py output')
    e.add_argument('--output', default='reports/isoform_rank_test/results')
    return p


if __name__ == '__main__':
    a = parser().parse_args()
    {'fetch': fetch, 'build': build, 'evaluate': evaluate}[a.command](a)
