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

Analysis extension: python -m unittest -v now runs **11 passing tests**.
Additional checks distinguish uniform scaling from reallocation, cover zero-mass
and disjoint-plan comparisons, type-specific transport/distance accounting,
matched references at fixed family/r0, changed inputs, output overwrite protection,
and distance-only runs without invented PPI comparisons.

analyze.py completed on the 45-case demo; its three figures were visually inspected.
Paper Methods and author-code notes are in NEXT_STEPS.md. The numerical pipeline
remains a proposed model, not a CCCvelo reproduction.

## Stage-one adapter — 2026-11-13

Final dependencies tested: Python 3.12, NumPy 2.3.5, SciPy 1.17.1,
pandas 2.2.3, Matplotlib 3.10.8, AnnData 0.12.19, h5py 3.15.1.
`pip check` passed. SciPy 1.17.0 was replaced because its sparse indexing
interface failed the backed-h5ad tests and is excluded by current AnnData.

17 tests passed. New coverage includes full-library normalization before LR
selection, preserving already processed expression, one-to-one symbol handling,
missing type coverage, unsupported complexes, pinned LR checksums, refusing
fractional data declared as counts, Delaunay exclusions/self-edge handling,
rejecting HTML downloads, and refusing changed prepared inputs.

The `fetch-lr` CLI actually downloaded and verified both CCCvelo CSV files from
commit fac213c600e7929b5d14e7d95dbdaabd722f9bc9 (3883 diffusion, 109 contact).
The inspect -> prepare -> run --sensitivity CLI completed nine OT settings and
45 distance cases using those real LR tables plus a **synthetic 16-spot h5ad
fixture**, with no PPI scores. The resulting distance figure was inspected.
These checks do not establish real tissue performance, an optimal r0, or
PLM-interact biological validity. The real mousebrain h5ad download and its
expression-layer provenance remain to be checked on the user's SSH server.

### Download URL correction — 2026-11-13

The user's download returned non-HDF5 content. Rechecking the scanned paper and
Spateo's original `mousebrain()` download dictionary established two transcription
errors in our default URL: `e5tu4drxda0lm0u` must be `c5tu4drxda01m0u`.
The replacement was matched programmatically to the literal in Spateo commit
1bd8a35e6e7b0dfc712d76b1c32987dc06652774, spateo/sample_data.py, line 157.
The cloud proxy denies Dropbox, so live success of the corrected data download
is not claimed; the source URL is verified and the SSH-server download remains next.

All 7 stage-one tests passed after this change, including a successful local-file
h5ad download and rejection/cleanup of a truncated HDF5 file that retains its
signature. The downloader now also checks Content-Length when present and readable
obs/var metadata before exposing the final file, and reports HTTP/content type on
invalid responses. No signature check was disabled to accept the HTML response.
