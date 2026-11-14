"""Summarize PLM-interact diagnostic scores for one human_pairs.py run (commot env).

Reads <run>/jobs.csv and every <run>/scores_<checkpoint>.csv whose .meta.json
records this run's plm_input.csv hash. Writes:

  pair_scores.csv   one row per pair x checkpoint: both orders and their mean
  group_summary.csv cognate vs cross per checkpoint and receptor-length group
  ligand_ranks.csv  for each cognate pair: rank of the cognate receptor among all
                    receptors in the set for that ligand (1 = scored highest)

Cross pairs are unannotated, not negatives; documented cross-reactive pairs are
reported but excluded from the AUROC.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from proteins import auroc, digest, write_json  # noqa: E402
from research import check  # noqa: E402

KEYS = ['ligand', 'receptor', 'set', 'receptor_group', 'ligand_len', 'receptor_len', 'known_cross_reactivity']


def load(run):
    run = Path(run)
    jobs = pd.read_csv(run/'jobs.csv')
    expected = digest(run/'plm_input.csv')
    frames = []
    for path in sorted(run.glob('scores_*.csv')):
        meta = json.loads(Path(str(path) + '.meta.json').read_text())
        check(meta['status'] == 'complete' and meta['input_sha256'] == expected,
              f'{path.name} was not scored from this run\'s plm_input.csv')
        check(digest(path) == meta['output_sha256'], f'{path.name} changed after scoring')
        s = pd.read_csv(path)
        check(set(s.pair_id) == set(jobs.pair_id), f'{path.name}: pair IDs differ from jobs.csv')
        m = jobs.merge(s[['pair_id', 'score']], on='pair_id', validate='one_to_one')
        lr = m[m.order == 'LR'][KEYS + ['score']].rename(columns={'score': 'score_LR'})
        rl = m[m.order == 'RL'][['ligand', 'receptor', 'score']].rename(columns={'score': 'score_RL'})
        wide = lr.merge(rl, on=['ligand', 'receptor'], validate='one_to_one')
        check(len(wide) * 2 == len(m), f'{path.name}: every pair needs exactly both orders')
        wide['score'] = (wide.score_LR + wide.score_RL) / 2
        wide.insert(0, 'checkpoint', meta['checkpoint_name'].split('/')[-1])
        wide['checkpoint_revision'] = meta['checkpoint_revision']
        frames.append(wide)
    check(frames, f'No scores_*.csv in {run}')
    return pd.concat(frames, ignore_index=True)


def summarize(scores):
    groups = []
    for (ckpt, group), d in scores.groupby(['checkpoint', 'receptor_group']):
        cog = d[d.set == 'cognate']
        cross = d[(d.set == 'cross') & d.known_cross_reactivity.isna()]
        groups.append(dict(checkpoint=ckpt, receptor_group=group, n_cognate=len(cog),
            cognate_median=cog.score.median(), cognate_min=cog.score.min(),
            cognate_frac_ge_0_5=(cog.score >= 0.5).mean() if len(cog) else None,
            n_cross=len(cross), cross_median=cross.score.median() if len(cross) else None,
            cross_frac_ge_0_5=(cross.score >= 0.5).mean() if len(cross) else None,
            auroc_cognate_vs_cross=auroc(cog.score, cross.score)))
    ranks = []
    for (ckpt, ligand), d in scores.groupby(['checkpoint', 'ligand']):
        d = d.sort_values('score', ascending=False).reset_index(drop=True)
        for i, row in d[d.set == 'cognate'].iterrows():
            ranks.append(dict(checkpoint=ckpt, ligand=ligand, receptor=row.receptor,
                              receptor_group=row.receptor_group, score=row.score,
                              rank_among_receptors=i + 1, n_receptors=len(d),
                              best_other_receptor=d.loc[d.index != i, 'receptor'].iloc[0] if len(d) > 1 else None))
    return pd.DataFrame(groups), pd.DataFrame(ranks)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run', required=True)
    args = p.parse_args()
    run = Path(args.run)
    scores = load(run)
    groups, ranks = summarize(scores)
    scores.to_csv(run/'pair_scores.csv', index=False)
    groups.to_csv(run/'group_summary.csv', index=False)
    ranks.to_csv(run/'ligand_ranks.csv', index=False)
    write_json(run/'summary.json', dict(status='complete',
        checkpoints=sorted(scores.checkpoint.unique()),
        inputs_sha256={f.name: digest(f) for f in sorted(run.glob('scores_*.csv'))},
        notes=['Mean of both input orders.', 'Cross pairs are unannotated, not negatives.',
               'Documented cross-reactive pairs are excluded from AUROC.',
               'Scores are model predictions, not measured affinities.']))
    pd.set_option('display.width', 220)
    cog = scores[scores.set == 'cognate'].pivot_table(index=['receptor_group', 'ligand', 'receptor'],
                                                      columns='checkpoint', values='score')
    print(cog.to_string(float_format=lambda v: f'{v:.3g}'))
    print(groups.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    if len(ranks):
        r = ranks.pivot_table(index=['ligand', 'receptor'], columns='checkpoint', values='rank_among_receptors')
        print(f'rank of cognate receptor per ligand (of {ranks.n_receptors.max()}):')
        print(r.to_string())


if __name__ == '__main__':
    main()
