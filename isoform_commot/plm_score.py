"""Score protein pairs with the official PLM-interact model (separate env).

Run in its own environment (Python 3.10, transformers==4.40.1, torch), see
requirements-plm.txt. This file deliberately imports nothing from the commot code.

  score     plm_input.csv (pair_id,query,text) -> scores.csv + scores.csv.meta.json
  selftest  sample of an official labelled test set (query,text,label) -> AUPR/AUROC,
            to check the installation and weights before trusting new scores

Differences from the official batch script, all deliberate:
  * one pair per forward pass, max_length = len(query)+len(text)+3 for THAT pair,
    and truncation is disabled and verified (the official script truncates silently);
  * output keyed by pair_id, with the checkpoint SHA256 recorded.

The PLMinteract class below is copied from the official repository
https://github.com/liudan111/PLM-interact (commit ab2ae6ae1aa81accf4c1f7f6f341164baf97809a,
PLMinteract/inference/inference_PPI_singleGPU.py), MIT License, Copyright (c) 2024 Dan LIU.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import pandas as pd

ESM_LETTERS = set('ACDEFGHIKLMNPQRSTVWYXBUZO')
OFFICIAL_MAX_TOKENS = 1603


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def load_jobs(path, columns=('pair_id', 'query', 'text')):
    jobs = pd.read_csv(path, dtype=str)
    missing = set(columns) - set(jobs)
    if missing:
        raise ValueError(f'{path} lacks columns {sorted(missing)}')
    if jobs[list(columns)].isna().any().any():
        raise ValueError('Empty sequence or pair_id')
    if 'pair_id' in columns and not jobs.pair_id.is_unique:
        raise ValueError('pair_id must be unique')
    for col in ['query', 'text']:
        jobs[col] = jobs[col].str.strip().str.upper()
        bad = set(''.join(jobs[col])) - ESM_LETTERS
        if bad:
            raise ValueError(f'Letters outside the ESM vocabulary in {col}: {sorted(bad)}')
    jobs['n_tokens'] = jobs['query'].str.len() + jobs['text'].str.len() + 3
    return jobs


def run_jobs(jobs, scorer, max_tokens=OFFICIAL_MAX_TOKENS, allow_long=False):
    too_long = jobs.loc[jobs.n_tokens > max_tokens, 'pair_id'].tolist()
    if too_long and not allow_long:
        raise ValueError(f'{len(too_long)} pairs exceed {max_tokens} tokens (official limit; training '
                         f'proteins were <=800 aa): {too_long[:10]}. Decide on a documented rule '
                         '(e.g. extracellular region) or pass --allow-long to score full length.')
    rows = []
    for i, job in enumerate(jobs.itertuples(index=False)):
        score, n_tokens = scorer(job.query, job.text)
        if n_tokens != job.n_tokens:
            raise ValueError(f'{job.pair_id}: tokenizer produced {n_tokens} tokens, expected '
                             f'{job.n_tokens}; refusing altered/truncated input')
        rows.append(dict(pair_id=getattr(job, 'pair_id', i), score=float(score), n_tokens=int(n_tokens)))
        if (i + 1) % 20 == 0:
            print(f'{i + 1}/{len(jobs)} pairs', flush=True)
    return pd.DataFrame(rows), too_long


def average_precision(labels, scores):
    order = np.argsort(-np.asarray(scores, float), kind='stable')
    y = np.asarray(labels, int)[order]
    if not y.sum():
        return None
    precision = np.cumsum(y) / np.arange(1, len(y) + 1)
    return float((precision * y).sum() / y.sum())


def auroc(labels, scores):
    y, s = np.asarray(labels, int), np.asarray(scores, float)
    pos, neg = s[y == 1], s[y == 0]
    if not len(pos) or not len(neg):
        return None
    return float(((pos[:, None] > neg[None, :]).sum() + 0.5*(pos[:, None] == neg[None, :]).sum())
                 / (len(pos)*len(neg)))


# ------------------------------------------------------------- model (torch only below)
def build_scorer(base_model, checkpoint, device=None):
    import torch
    import torch.nn as nn
    import torch.nn.functional as F
    import transformers
    from transformers import AutoConfig, AutoModelForMaskedLM, AutoTokenizer

    class PLMinteract(nn.Module):  # official code, see module docstring
        def __init__(self, model_name, num_labels, config, device, embedding_size):
            super(PLMinteract, self).__init__()
            self.esm_mask = AutoModelForMaskedLM.from_pretrained(model_name, config=config)
            self.embedding_size = embedding_size
            self.classifier = nn.Linear(embedding_size, 1)  # embedding_size
            self.num_labels = num_labels
            self.device = device

        def forward_test(self, features):
            embedding_output = self.esm_mask.base_model(**features, return_dict=True)
            embedding = embedding_output.last_hidden_state[:, 0, :]  # cls token
            embedding = F.relu(embedding)
            logits = self.classifier(embedding)
            logits = logits.view(-1)
            probability = torch.sigmoid(logits)
            return probability

    device = torch.device(device or ('cuda:0' if torch.cuda.is_available() else 'cpu'))
    config = AutoConfig.from_pretrained(base_model)
    tokenizer = AutoTokenizer.from_pretrained(base_model)
    state = torch.load(checkpoint, map_location='cpu', weights_only=True)
    embedding_size = int(state['classifier.weight'].shape[1])
    if embedding_size != config.hidden_size:
        raise ValueError(f'Checkpoint head ({embedding_size}) does not match base model '
                         f'hidden size ({config.hidden_size}); use the matching ESM-2 base model')
    model = PLMinteract(base_model, 1, config=config, device=device, embedding_size=embedding_size)
    model.load_state_dict(state, strict=True)
    model.eval().to(device)

    def scorer(query, text):
        features = tokenizer(query, text, padding=False, truncation=False, return_tensors='pt')
        n_tokens = int(features['input_ids'].shape[1])
        features = {k: v.to(device) for k, v in features.items()}
        with torch.no_grad():
            probability = model.forward_test(features)
        return float(probability.item()), n_tokens

    info = dict(torch=torch.__version__, transformers=transformers.__version__,
                device=str(device), embedding_size=embedding_size,
                cuda=torch.version.cuda if torch.cuda.is_available() else None)
    return scorer, info


def model_meta(args, info):
    return dict(checkpoint=str(Path(args.checkpoint).resolve()), checkpoint_sha256=digest(args.checkpoint),
                checkpoint_name=args.checkpoint_name, checkpoint_revision=args.checkpoint_revision,
                base_model=args.base_model, python=platform.python_version(), **info)


def score(args):
    out = Path(args.output)
    meta_path = Path(str(out) + '.meta.json')
    if out.exists() or meta_path.exists():
        raise FileExistsError(f'{out} exists; use a new output name')
    jobs = load_jobs(args.jobs)
    scorer, info = build_scorer(args.base_model, args.checkpoint, args.device)
    started = now()
    result, too_long = run_jobs(jobs, scorer, args.max_tokens, args.allow_long)
    result.to_csv(out, index=False)
    meta_path.write_text(json.dumps(dict(status='complete', input=str(Path(args.jobs).resolve()),
        input_sha256=digest(args.jobs), n_pairs=len(result), started_utc=started, finished_utc=now(),
        truncation='disabled; token count verified per pair', max_tokens=args.max_tokens,
        scored_beyond_max_tokens=too_long, output_sha256=digest(out), **model_meta(args, info)), indent=2))
    print(f'Scored {len(result)} pairs -> {out}')


def selftest(args):
    out = Path(args.output)
    if out.exists():
        raise FileExistsError(f'{out} exists')
    data = pd.read_csv(args.test_file, dtype={'query': str, 'text': str})
    if 'label' not in data:
        raise ValueError('Self-test file needs query,text,label')
    data = data.sample(n=min(args.n, len(data)), random_state=args.seed).reset_index(drop=True)
    data['pair_id'] = [f'selftest_{i}' for i in range(len(data))]
    path = out.with_suffix('.input.csv')
    data[['pair_id', 'query', 'text']].to_csv(path, index=False)
    jobs = load_jobs(path)
    scorer, info = build_scorer(args.base_model, args.checkpoint, args.device)
    result, _ = run_jobs(jobs, scorer, args.max_tokens, allow_long=False)
    labels = data.label.to_numpy(int)
    report = dict(status='complete', test_file=str(Path(args.test_file).resolve()),
        test_file_sha256=digest(args.test_file), n=len(result), positives=int(labels.sum()), seed=args.seed,
        aupr=average_precision(labels, result.score), auroc=auroc(labels, result.score),
        note='Compare with the cross-species result reported in the PLM-interact paper for the same species; '
             'a large gap suggests a wrong base model, checkpoint or environment.', **model_meta(args, info))
    out.write_text(json.dumps(report, indent=2))
    print(json.dumps({k: report[k] for k in ['n', 'positives', 'aupr', 'auroc']}, indent=2))


def parser():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    for name in ['score', 'selftest']:
        s = sub.add_parser(name)
        s.add_argument('--checkpoint', required=True, help='PLM-interact pytorch_model.bin')
        s.add_argument('--base-model', required=True,
                       help='local dir or HF name of the matching ESM-2, e.g. facebook/esm2_t33_650M_UR50D')
        s.add_argument('--checkpoint-name', required=True, help='e.g. danliu1226/PLM-interact-650M-humanV11')
        s.add_argument('--checkpoint-revision', required=True, help='HF commit hash of the checkpoint')
        s.add_argument('--device'); s.add_argument('--output', required=True)
        s.add_argument('--max-tokens', type=int, default=OFFICIAL_MAX_TOKENS)
    sub.choices['score'].add_argument('--jobs', required=True, help='plm_input.csv from proteins.py pairs')
    sub.choices['score'].add_argument('--allow-long', action='store_true')
    t = sub.choices['selftest']
    t.add_argument('--test-file', required=True, help='e.g. data/test/mouse.ppi.qrels.seq.test.csv')
    t.add_argument('--n', type=int, default=2000); t.add_argument('--seed', type=int, default=2)
    return p


if __name__ == '__main__':
    arguments = parser().parse_args()
    score(arguments) if arguments.command == 'score' else selftest(arguments)
