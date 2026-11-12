# Validation

2026-11-12: Python 3.12.14, NumPy 2.3.5, SciPy 1.17.0,
pandas 2.2.3, Matplotlib 3.10.8.

`python -m unittest -v`: **6 tests passed**.

- Coordinate-unit invariance and the small-r0 limit of the reference-normalized kernel.
- Valid zero-cost contact edges, spatial exclusions, and rejected invalid scores.
- Zero transport on completely disconnected candidates.
- Changed allocation under a controlled PPI cost perturbation with competition.
- Official L1-UNOT solution compared to independent constrained primal optimization (SLSQP).
- Full synthetic sweep, refusal to overwrite inputs/results, and rejected missing PPI coverage.

These are software/numerical checks, not biological validation.
The demo uses synthetic expression and invented scores, not PLM-interact predictions.
CCCvelo-specific preprocessing, real-data performance and optimal r0 remain unverified.
