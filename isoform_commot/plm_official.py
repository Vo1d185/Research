"""Run the OFFICIAL PLM-interact inference script unchanged and key its output by pair_id.

Runs in the plminteract env (prepare, collect) and launches the official script with the
separate venv that adds the official imports (sentence-transformers, datasets, scikit-learn),
see reports/plm_lr_benchmark_plan.md section 3.3.

  prepare   jobs CSV (pair_id,query,text) -> <dir>/official_input.csv (query,text: the official
            format) and <dir>/official_ids.csv (same row order) with the token count of each pair
            before and after the official truncation, obtained with the same tokenizer call as the
            official collate (truncation='longest_first', max_length)
  run       python -m torch.distributed.run --standalone --nproc_per_node=1 inference_PPI_singleGPU.py
            with the README single-GPU parameters; TORCH_FORCE_WEIGHTS_ONLY_LOAD=1 makes the
            script's torch.load weights-only (the script itself is not modified)
  collect   <dir>/pred_scores.csv (row order) -> scores CSV keyed by pair_id + .meta.json
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess

import pandas as pd

OFFICIAL_COMMIT = 'ab2ae6ae1aa81accf4c1f7f6f341164baf97809a'
README_ARGS = dict(seed='2', batch_size_val='16', max_length='1603', embedding_size='1280',
                   model_name='esm2_t33_650M_UR50D')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def prepare(args):
    from transformers import AutoTokenizer
    out = Path(args.dir)
    out.mkdir(parents=True)
    jobs = pd.read_csv(args.jobs, dtype=str)
    if not {'pair_id', 'query', 'text'} <= set(jobs) or not jobs.pair_id.is_unique:
        raise ValueError('jobs need unique pair_id, query, text')
    for col in ['query', 'text']:
        jobs[col] = jobs[col].str.strip().str.upper()
    tok = AutoTokenizer.from_pretrained(args.base_model)
    max_length = int(README_ARGS['max_length'])
    rows = []
    for job in jobs.itertuples(index=False):
        full = tok(job.query, job.text, padding=False, truncation=False)['input_ids']
        used = tok(job.query, job.text, padding=False, truncation='longest_first', max_length=max_length)
        sep = used['input_ids'].index(tok.eos_token_id)
        rows.append(dict(pair_id=job.pair_id, query_len=len(job.query), text_len=len(job.text),
                         n_tokens_full=len(full), n_tokens_used=len(used['input_ids']),
                         query_len_used=sep - 1, text_len_used=len(used['input_ids']) - sep - 2,
                         truncated=len(full) > max_length))
    ids = pd.DataFrame(rows)
    if (ids.n_tokens_full != ids.query_len + ids.text_len + 3).any():
        raise ValueError('Unexpected tokenisation (expected one token per residue + 3 special tokens)')
    jobs[['query', 'text']].to_csv(out/'official_input.csv', index=False)
    ids.to_csv(out/'official_ids.csv', index=False)
    (out/'prepare.json').write_text(json.dumps(dict(jobs=str(Path(args.jobs).resolve()),
        jobs_sha256=digest(args.jobs), n=len(ids), truncated=int(ids.truncated.sum()),
        official_input_sha256=digest(out/'official_input.csv'), max_length=max_length,
        base_model=args.base_model, created_utc=now()), indent=2))
    print(f'{len(ids)} pairs, {int(ids.truncated.sum())} truncated by the official max_length -> {out}')


def run(args):
    out = Path(args.dir).resolve()
    src = Path(args.official_src).resolve()
    head = subprocess.run(['git', 'rev-parse', 'HEAD'], cwd=src, capture_output=True, text=True).stdout.strip()
    if head != OFFICIAL_COMMIT:
        raise ValueError(f'official source at {head}, expected {OFFICIAL_COMMIT}')
    changed = [line for line in subprocess.run(['git', 'status', '--porcelain', '--untracked-files=no'], cwd=src,
               capture_output=True, text=True).stdout.splitlines() if '__pycache__/' not in line]
    if changed:  # the repository tracks bytecode caches; those are not source changes
        raise ValueError(f'official source has local modifications: {changed}')
    if (out/'pred_scores.csv').exists():
        raise FileExistsError(out/'pred_scores.csv')
    base = Path(args.base_model).resolve()
    command = [args.venv_python, '-m', 'torch.distributed.run', '--standalone', '--nproc_per_node=1',
               str(src/'PLMinteract'/'inference'/'inference_PPI_singleGPU.py'),
               '--seed', README_ARGS['seed'], '--batch_size_val', README_ARGS['batch_size_val'],
               '--test_filepath', str(out/'official_input.csv'), '--model_name', README_ARGS['model_name'],
               '--embedding_size', README_ARGS['embedding_size'], '--output_filepath', str(out) + '/',
               '--resume_from_checkpoint', str(Path(args.checkpoint).resolve()),
               '--max_length', README_ARGS['max_length'], '--offline_model_path', str(base.parent) + '/']
    if base.name != README_ARGS['model_name']:
        raise ValueError(f'base model directory must be named {README_ARGS["model_name"]}')
    env = dict(os.environ, PYTHONPATH=str(src), TORCH_FORCE_WEIGHTS_ONLY_LOAD='1', PYTHONNOUSERSITE='1',
               PYTHONDONTWRITEBYTECODE='1')
    started = now()
    subprocess.run(command, cwd=out, env=env, check=True)
    (out/'run.json').write_text(json.dumps(dict(command=command, official_commit=head, started_utc=started,
        finished_utc=now(), checkpoint=str(Path(args.checkpoint).resolve()), checkpoint_sha256=digest(args.checkpoint),
        checkpoint_name=args.checkpoint_name, checkpoint_revision=args.checkpoint_revision,
        env=dict(TORCH_FORCE_WEIGHTS_ONLY_LOAD='1', PYTHONPATH=str(src)),
        note='official script unmodified; launcher flag --standalone only picks a free rendezvous port'), indent=2))
    print(f'official scores -> {out}/pred_scores.csv')


def collect(args):
    out = Path(args.dir)
    ids = pd.read_csv(out/'official_ids.csv', dtype={'pair_id': str})
    pred = pd.read_csv(out/'pred_scores.csv', header=None)
    if pred.shape != (len(ids), 1):
        raise ValueError(f'pred_scores.csv has shape {pred.shape}, expected ({len(ids)}, 1)')
    ids['score'] = pred[0].astype(float).to_numpy()
    result = Path(args.output)
    if result.exists():
        raise FileExistsError(result)
    ids.to_csv(result, index=False)
    run_meta = json.loads((out/'run.json').read_text())
    Path(str(result) + '.meta.json').write_text(json.dumps(dict(status='complete', scorer='official',
        prepare=json.loads((out/'prepare.json').read_text()), run=run_meta,
        pred_scores_sha256=digest(out/'pred_scores.csv'), output_sha256=digest(result)), indent=2))
    print(f'{len(ids)} scores ({int(ids.truncated.sum())} truncated) -> {result}')


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    a = sub.add_parser('prepare')
    a.add_argument('--jobs', required=True); a.add_argument('--dir', required=True)
    a.add_argument('--base-model', default='models/esm2_t33_650M_UR50D')
    b = sub.add_parser('run')
    b.add_argument('--dir', required=True); b.add_argument('--checkpoint', required=True)
    b.add_argument('--checkpoint-name', required=True); b.add_argument('--checkpoint-revision', required=True)
    b.add_argument('--base-model', default='models/esm2_t33_650M_UR50D')
    b.add_argument('--official-src', default='models/PLM-interact-src')
    b.add_argument('--venv-python', default='/home/liu/venvs/plminteract-official/bin/python')
    c = sub.add_parser('collect')
    c.add_argument('--dir', required=True); c.add_argument('--output', required=True)
    return p


if __name__ == '__main__':
    arguments = parser().parse_args()
    {'prepare': prepare, 'run': run, 'collect': collect}[arguments.command](arguments)
