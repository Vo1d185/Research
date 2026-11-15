"""Random 2D fields for the PDE r0-lambda benchmark (reports/pde_r0_plan.md, section 4.2) and their statistics.

A field is f = min(sum_k G_k, 1) with values below 0.1 set to 0, where each G_k is an anisotropic Gaussian of
peak 1 (centre uniform on the domain, both axis widths uniform on [sigma_lo, sigma_hi], angle uniform on
[0, 180) degrees) and the number of blobs K is uniform on {1, ..., k_max}. Draw order per field: K, then per
blob: centre x, centre y, width 1, width 2, angle. A scenario (one binding structure) draws its ligand
production fields first, then its receptor initial fields, from numpy.random.default_rng([structure, repeat]).

Statistics used for calibration (computed identically on archived and generated fields):
  support   fraction of grid points with f > 0
  blobs     number of 4-connected components of {f > 0} (scipy.ndimage.label default)
  overlap   per interacting pair (i, j): |S(p_i) & S(R0_j)| / |S(R0_j)|
"""
from dataclasses import asdict, dataclass

import numpy as np
from scipy import ndimage
from scipy.stats import ks_2samp

from simulation.pde import gaussian2d

THRESHOLD = 0.1
DOMAIN = (-200.0, 200.0)
CALIBRATION_GRID = dict(k_max=list(range(1, 11)), sigma_hi=[20.0, 30.0, 40.0, 50.0, 60.0, 80.0],
                        sigma_ratio=[0.25, 0.5])
CALIBRATION_REPEATS = range(1001, 1011)       # seeds [structure, 1000 + s], s = 1..10; disjoint from 1..5
STATISTICS = ('support', 'blobs', 'overlap')


@dataclass(frozen=True)
class FieldParams:
    k_max: int
    sigma_hi: float
    sigma_ratio: float            # sigma_lo / sigma_hi

    @property
    def sigma_lo(self):
        return self.sigma_ratio * self.sigma_hi


def grid_points(n=50, lo=DOMAIN[0], hi=DOMAIN[1]):
    """Cell centres (n, n, 2), as the archived grid.npy (simulation.pde.Grid)."""
    c = lo + (hi - lo) / n * (np.arange(n) + 0.5)
    return np.stack(np.meshgrid(c, c, indexing='ij'), axis=-1)


def draw_field(rng, xy, params):
    k = int(rng.integers(1, params.k_max + 1))
    total = np.zeros(xy.shape[:-1])
    for _ in range(k):
        cx, cy = rng.uniform(*DOMAIN), rng.uniform(*DOMAIN)
        s1, s2 = rng.uniform(params.sigma_lo, params.sigma_hi), rng.uniform(params.sigma_lo, params.sigma_hi)
        angle = rng.uniform(0.0, 180.0)
        total += gaussian2d(xy, (cx, cy), (s1, s2), angle)
    f = np.minimum(total, 1.0)
    f[f < THRESHOLD] = 0.0
    return f


def scenario_fields(structure, repeat, n_ligands, n_receptors, params, xy=None):
    """(production (n_l, n, n), receptor0 (n_r, n, n)) from default_rng([structure, repeat])."""
    xy = grid_points() if xy is None else xy
    rng = np.random.default_rng([int(structure), int(repeat)])
    p = np.stack([draw_field(rng, xy, params) for _ in range(n_ligands)])
    r = np.stack([draw_field(rng, xy, params) for _ in range(n_receptors)])
    return p, r


def support(f):
    return float(np.mean(f > 0))


def blobs(f):
    return int(ndimage.label(f > 0)[1])


def overlap(p, r):
    s = r > 0
    return float((s & (p > 0)).sum() / s.sum()) if s.any() else np.nan


def scenario_statistics(p, r, binding):
    """Per-field support and blob counts, per-pair overlaps (binding: n_l x n_r, > 0 = interacting)."""
    fields = list(p) + list(r)
    pairs = [(i, j) for i in range(binding.shape[0]) for j in range(binding.shape[1]) if binding[i, j] > 0]
    return dict(support=[support(f) for f in fields], blobs=[blobs(f) for f in fields],
                overlap=[overlap(p[i], r[j]) for i, j in pairs])


def pooled(stats_list):
    return {k: np.array([v for s in stats_list for v in s[k]], float) for k in STATISTICS}


def ks_objective(reference, sample):
    """Sum over the three statistics of the two-sample KS statistic, and the parts."""
    parts = {k: float(ks_2samp(reference[k], sample[k]).statistic) for k in STATISTICS}
    return sum(parts.values()), parts


def calibration_grid():
    """All parameter combinations in grid order (k_max, then sigma_hi, then ratio, each ascending)."""
    g = CALIBRATION_GRID
    return [FieldParams(k, s, r) for k in g['k_max'] for s in g['sigma_hi'] for r in g['sigma_ratio']]


def calibrate(reference, structures, xy=None):
    """reference: pooled archived statistics; structures: {structure: binding matrix}.
    Returns (rows for every combination, chosen FieldParams). Ties: first in grid order."""
    xy = grid_points() if xy is None else xy
    rows = []
    for params in calibration_grid():
        stats = [scenario_statistics(*scenario_fields(c, s, a.shape[0], a.shape[1], params, xy), a)
                 for c, a in sorted(structures.items()) for s in CALIBRATION_REPEATS]
        total, parts = ks_objective(reference, pooled(stats))
        rows.append(dict(asdict(params), objective=total, **{f'ks_{k}': v for k, v in parts.items()}))
    best = min(range(len(rows)), key=lambda i: (rows[i]['objective'], i))
    return rows, calibration_grid()[best]
