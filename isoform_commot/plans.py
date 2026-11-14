"""Transport-plan files: dense (stage 1) and sparse COO (stage 2, COMMOT baseline).

Sparse case_*.npz keys
  format='coo', shape=(n_lig*n, n_rec*n), row, col  (eligible entries only)
  C (cost or eligibility record on those entries), T (transport on those entries)
  a, b, spots, ligands, receptors
Row index = ligand_index*n + sender_spot; column index = receptor_index*n + receiver_spot.
Entries that are not stored are forbidden (infinite cost) and carry no transport.
"""
from pathlib import Path

import numpy as np
from scipy import sparse

from research import check


def save_sparse(path, row, col, C, T, shape, a, b, spots, ligands, receptors, **extra):
    row, col = np.asarray(row, np.int64), np.asarray(col, np.int64)
    C, T = np.asarray(C, float), np.asarray(T, float)
    check(row.shape == col.shape == C.shape == T.shape, 'Entry arrays differ in length')
    check(len(shape) == 2 and (row < shape[0]).all() and (col < shape[1]).all(), 'Entry outside matrix')
    check(np.isfinite(C).all() and np.isfinite(T).all() and (T >= 0).all(), 'Invalid stored cost/transport')
    np.savez_compressed(path, format='coo', shape=np.asarray(shape, np.int64), row=row, col=col, C=C, T=T,
                        a=np.asarray(a, float), b=np.asarray(b, float), spots=np.asarray(spots, str),
                        ligands=np.asarray(ligands, str), receptors=np.asarray(receptors, str), **extra)


def load(path):
    """Return dict with 'sparse' flag; dense files keep T/C as arrays, sparse ones as entry arrays."""
    with np.load(Path(path), allow_pickle=False) as z:
        out = {k: z[k] for k in z.files}
    out['sparse'] = 'format' in out and str(out['format']) == 'coo'
    if out['sparse']:
        out['shape'] = tuple(int(v) for v in out['shape'])
    return out


def to_coo(plan, key='T'):
    if plan['sparse']:
        return sparse.coo_matrix((plan[key], (plan['row'], plan['col'])), shape=plan['shape'])
    M = np.where(np.isfinite(plan['C']), plan[key], 0.0) if key == 'C' else plan[key]
    return sparse.coo_matrix(M)


def received_signal(plan, pairs, n):
    """Per LR pair and per receptor: received transport per receiver spot (sum over senders)."""
    ligands, receptors = list(plan['ligands'].astype(str)), list(plan['receptors'].astype(str))
    li, ri = {v: i for i, v in enumerate(ligands)}, {v: i for i, v in enumerate(receptors)}
    keys = [(row.ligand, row.receptor) for row in pairs.itertuples(index=False)]
    if plan['sparse']:
        pair_of = np.full((len(ligands), len(receptors)), -1, np.int64)
        for p, (l, r) in enumerate(keys):
            pair_of[li[l], ri[r]] = p
        p = pair_of[plan['row'] // n, plan['col'] // n]
        check((plan['T'][p < 0] <= 1e-12).all(), 'Transport on a block that is not an LR pair')
        keep = p >= 0
        flat = np.bincount(p[keep] * n + plan['col'][keep] % n, weights=plan['T'][keep],
                           minlength=len(keys) * n).reshape(len(keys), n)
        per_pair = {k: flat[i] for i, k in enumerate(keys)}
    else:
        T = plan['T']
        per_pair = {(l, r): T[li[l]*n:(li[l]+1)*n, ri[r]*n:(ri[r]+1)*n].sum(axis=0) for l, r in keys}
    per_receptor = {}
    for (ligand, receptor), signal in per_pair.items():
        per_receptor[receptor] = per_receptor.get(receptor, 0) + signal
    return per_pair, per_receptor
