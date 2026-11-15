"""Reaction-diffusion model of ligand-receptor binding (COMMOT Supplementary Note 2.1, Eq. 19; Extended Data Fig. 1a).

    d[L_i]/dt    = D_i lap[L_i] - sum_j a_ij [L_i][R_j] + sum_j b_ij [L_iR_j] - c_i [L_i] + p_i(x)
    d[L_iR_j]/dt = a_ij [L_i][R_j] - b_ij [L_iR_j]          (- e_ij [L_iR_j], archive field, 0 there)
    d[R_j]/dt    = - sum_i a_ij [L_i][R_j] + sum_i b_ij [L_iR_j]   (- g_j [R_j], archive field, 0 there)

p_i(x) is a constant ligand production field. Eq. 19 and the figure omit it, but COMMOT's archived simulations
start from [L_i] = 0 and accumulate ligand from the files ligand_prod_i.npy (reports/pde_reimplementation.md).
Receptors start from receptor_init; complexes start at 0. No-flux (Neumann) boundaries on a cell-centred grid,
as in py-pde's CartesianGrid; receptors and complexes do not diffuse.

Steady state (e = g = 0): binding balances (a L R = b C), so the free ligand solves the linear problem
(c_i - D_i lap) L_i = p_i, and with K_ij = a_ij / b_ij the receptor is shared as
R_j = R0_j / (1 + sum_i K_ij L_i), C_ij = K_ij L_i R_j. steady_state() uses this; integrate() is independent.
"""
from dataclasses import dataclass, field

import numpy as np
from scipy import sparse
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import spsolve


@dataclass
class Grid:
    """Cell-centred grid on [lo, hi] (1D) or [lo, hi]^2 (2D) with n cells per axis."""
    n: int
    lo: float
    hi: float
    dim: int = 1

    @property
    def dx(self):
        return (self.hi - self.lo) / self.n

    @property
    def centres(self):
        return self.lo + self.dx * (np.arange(self.n) + 0.5)

    @property
    def shape(self):
        return (self.n,) * self.dim

    def coords(self):
        c = self.centres
        return c if self.dim == 1 else np.stack(np.meshgrid(c, c, indexing='ij'), axis=-1)

    def laplacian(self, u):
        """Second differences with zero-flux ghost cells (ghost = edge value)."""
        p = np.pad(u, 1, mode='edge')
        if self.dim == 1:
            return (p[2:] + p[:-2] - 2*u) / self.dx**2
        return (p[2:, 1:-1] + p[:-2, 1:-1] + p[1:-1, 2:] + p[1:-1, :-2] - 4*u) / self.dx**2

    def laplacian_matrix(self):
        one = sparse.diags([np.ones(self.n - 1), -2*np.ones(self.n), np.ones(self.n - 1)], [-1, 0, 1], format='lil')
        one[0, 0] = one[-1, -1] = -1                      # zero flux
        one = one.tocsr() / self.dx**2
        if self.dim == 1:
            return one
        eye = sparse.identity(self.n, format='csr')
        return (sparse.kron(one, eye) + sparse.kron(eye, one)).tocsr()


@dataclass
class Model:
    """Rates as in COMMOT's archived meta_data: D (n_l), a, b (n_l x n_r; a = 0: no binding), c (n_l)."""
    D: np.ndarray
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    e: np.ndarray = None          # complex degradation (archive field; 0)
    g: np.ndarray = None          # receptor degradation (archive field; 0)
    pairs: list = field(init=False)

    def __post_init__(self):
        self.D, self.a, self.b, self.c = (np.asarray(v, float) for v in (self.D, self.a, self.b, self.c))
        nl, nr = self.a.shape
        self.e = np.zeros((nl, nr)) if self.e is None else np.asarray(self.e, float)
        self.g = np.zeros(nr) if self.g is None else np.asarray(self.g, float)
        assert self.D.shape == self.c.shape == (nl,) and self.b.shape == self.e.shape == (nl, nr)
        self.pairs = [(i, j) for i in range(nl) for j in range(nr) if self.a[i, j] > 0]

    @property
    def nl(self):
        return self.a.shape[0]

    @property
    def nr(self):
        return self.a.shape[1]


class System:
    """State y = [L_1..L_nl, R_1..R_nr, C_pair for pair in model.pairs], each a field on the grid."""

    def __init__(self, model, grid, production, receptor0):
        self.m, self.grid = model, grid
        self.p = np.asarray(production, float).reshape((model.nl,) + grid.shape)
        self.r0 = np.asarray(receptor0, float).reshape((model.nr,) + grid.shape)
        self.size = int(np.prod(grid.shape))
        self.nfields = model.nl + model.nr + len(model.pairs)

    def split(self, y):
        f = y.reshape((self.nfields,) + self.grid.shape)
        nl, nr = self.m.nl, self.m.nr
        return f[:nl], f[nl:nl + nr], f[nl + nr:]

    def initial(self):
        y = np.zeros((self.nfields,) + self.grid.shape)
        y[self.m.nl:self.m.nl + self.m.nr] = self.r0
        return y.ravel()

    def rhs(self, t, y):
        L, R, C = self.split(y)
        dL = np.stack([self.m.D[i]*self.grid.laplacian(L[i]) - self.m.c[i]*L[i] + self.p[i] for i in range(self.m.nl)])
        dR = -self.m.g.reshape((-1,) + (1,)*self.grid.dim) * R
        dC = np.empty_like(C)
        for k, (i, j) in enumerate(self.m.pairs):
            flux = self.m.a[i, j]*L[i]*R[j] - self.m.b[i, j]*C[k]
            dL[i] -= flux
            dR[j] -= flux
            dC[k] = flux - self.m.e[i, j]*C[k]
        return np.concatenate([dL.ravel(), dR.ravel(), dC.ravel()])

    def integrate(self, times, rtol=1e-8, atol=1e-10, method='RK45', y0=None):
        """Fields at the requested times; method/rtol/atol passed to scipy.integrate.solve_ivp."""
        sol = solve_ivp(self.rhs, (times[0], times[-1]), self.initial() if y0 is None else y0, t_eval=times,
                        method=method, rtol=rtol, atol=atol)
        if not sol.success:
            raise RuntimeError(sol.message)
        return [self.split(sol.y[:, k]) for k in range(len(times))], sol

    def steady_state(self):
        """Closed form for e = g = 0 (module docstring)."""
        assert not self.m.e.any() and not self.m.g.any(), 'closed form needs no complex/receptor degradation'
        lap = self.grid.laplacian_matrix()
        eye = sparse.identity(self.size, format='csr')
        L = np.stack([spsolve((self.m.c[i]*eye - self.m.D[i]*lap).tocsc(), self.p[i].ravel()).reshape(self.grid.shape)
                      for i in range(self.m.nl)])
        K = np.divide(self.m.a, self.m.b, out=np.zeros_like(self.m.a), where=self.m.a > 0)
        denom = 1 + np.einsum('ij,i...->j...', K, L)
        R = self.r0 / denom
        C = np.stack([K[i, j]*L[i]*R[j] for i, j in self.m.pairs])
        return L, R, C

    def run_to_steady_state(self, t_chunk=500.0, tol=1e-8, t_max=1e5, method='LSODA', rtol=1e-10, atol=1e-12):
        """Integrate in chunks until the largest change of any complex over a chunk, relative to the largest
        complex value, is below tol. Independent of steady_state()."""
        y, t, history = self.initial(), 0.0, []
        while t < t_max:
            sol = solve_ivp(self.rhs, (t, t + t_chunk), y, method=method, rtol=rtol, atol=atol)
            if not sol.success:
                raise RuntimeError(sol.message)
            y_new, t = sol.y[:, -1], t + t_chunk
            C_old, C_new = self.split(y)[2], self.split(y_new)[2]
            change = float(np.abs(C_new - C_old).max() / max(np.abs(C_new).max(), 1e-300))
            history.append((t, change))
            y = y_new
            if change < tol:
                return self.split(y), t, history
        raise RuntimeError(f'no steady state by t = {t_max} (last relative change {change:.2e})')


def bump(x, centre, width):
    """Raised-cosine bump of half-width `width` (smooth, compact support), maximum 1."""
    r = np.abs(np.asarray(x, float) - centre) / width
    return np.where(r < 1, 0.5*(1 + np.cos(np.pi*r)), 0.0)


def gaussian2d(xy, centre, sigma, angle_deg=0.0):
    """Anisotropic 2D Gaussian (maximum 1); sigma = (s_major, s_minor), rotated by angle_deg."""
    t = np.deg2rad(angle_deg)
    d = xy - np.asarray(centre, float)
    u = d[..., 0]*np.cos(t) + d[..., 1]*np.sin(t)
    v = -d[..., 0]*np.sin(t) + d[..., 1]*np.cos(t)
    return np.exp(-0.5*((u/sigma[0])**2 + (v/sigma[1])**2))


def crescent(xy, centre, radius, width, start_deg, end_deg):
    """Arc of an annulus (radial Gaussian profile) between two polar angles (degrees, counter-clockwise)."""
    d = xy - np.asarray(centre, float)
    r = np.hypot(d[..., 0], d[..., 1])
    ang = np.rad2deg(np.arctan2(d[..., 1], d[..., 0])) % 360
    span = (end_deg - start_deg) % 360
    inside = ((ang - start_deg) % 360) <= span
    taper = np.sin(np.pi*np.clip(((ang - start_deg) % 360)/span, 0, 1))**0.5
    return np.exp(-0.5*((r - radius)/width)**2) * inside * taper
