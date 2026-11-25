"""
Sparse PDL VAR/VECM for high-resolution price data (Python port of MVARi.m).

Setup (as in Hasbrouck's MVARi with polynomial distributed lags of degree 0):
    p_j(t)   price level j on a grid of 1/second_factor seconds (step function)
    dp_j(t)  = p_j(t) - p_j(t-1), nonzero only when p_j changes
    X_jk(t)  = sum_{l=a_k+1}^{a_k+n_k} dp_j(t-l)     (PDL term k; a_k = k_offset)
    e_m(t)   = p_1(t) - p_{m+1}(t)                    (error-correction terms)

    dp_i(t) = sum_{j,k} beta_{ijk} X_jk(t) + sum_m gamma_im (e_m(t-1) - mean_m) + u_i(t),
    estimated by OLS for t in [first, last].

At 10 microsecond resolution a trading day has ~2.25 billion grid points, so the
data matrix z is never formed. Hasbrouck builds z'z from sparse cross-product
routines (Crossproduct.m). Here we use an equivalent exact shortcut: every
column of z is piecewise constant, changing only at "breakpoints" (price-change
times shifted by 0, 1 and each PDL window edge). Between consecutive breakpoints
z(t) is constant, so

    z'z = sum_s  w_s * z(b_s) z(b_s)'

where b_s are the breakpoints and w_s the number of grid points until the next
one. The dp columns are nonzero only at a change time tau, and tau+1 is always a
breakpoint, so those segments have length one.
"""
import math
from dataclasses import dataclass
from typing import List, Optional, Sequence

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Poly:
    """Degree-zero polynomial distributed lag: sums dp(t-l) for l = k_offset+1 ... k_offset+n."""
    name: str
    n: int
    k_offset: int

    @property
    def var_name(self) -> str:
        return f"{self.name}d0"


def pdl_polys(second_factor: int, lag_pdl_max: int = 10 ** 6) -> List[Poly]:
    """The PDL specification used in MVARiPartSIP01.m (and the other 01/02 programs):
    p0 = lag 1, p1 = lags 2-10, p2 = lags 11-100, ... up to 10 seconds of lags."""
    polys = [Poly("p0", 1, 0)]
    top = min(math.floor(math.log10(second_factor) + 1e-9),
              math.floor(math.log10(lag_pdl_max) + 1e-9) - 1)
    for i in range(top + 1):
        polys.append(Poly(f"p{i + 1}", 9 * 10 ** i, 10 ** i))
    return polys


@dataclass
class PriceSeries:
    """A price level that persists until changed (Hasbrouck's spl class).
    t: strictly increasing grid indices of the changes; v: the new price at each."""
    name: str
    t: np.ndarray
    v: np.ndarray

    def __post_init__(self):
        self.t = np.asarray(self.t, dtype=np.int64)
        self.v = np.asarray(self.v, dtype=np.float64)
        if len(self.t) != len(self.v):
            raise ValueError(f"{self.name}: t and v differ in length")
        if len(self.t) < 2:
            raise ValueError(f"{self.name}: fewer than two observations")
        if np.any(np.diff(self.t) <= 0):
            raise ValueError(f"{self.name}: t must be strictly increasing")

    def level(self, q: np.ndarray) -> np.ndarray:
        """Price in effect at grid times q (NaN before the first observation)."""
        idx = np.searchsorted(self.t, q, side="right") - 1
        out = np.full(q.shape, np.nan)
        ok = idx >= 0
        out[ok] = self.v[idx[ok]]
        return out

    def cum_dp(self, q: np.ndarray) -> np.ndarray:
        """Cumulative price change up to q (zero before the first observation).
        Sums of dp over a window are differences of this function."""
        idx = np.searchsorted(self.t, q, side="right") - 1
        return np.where(idx >= 0, self.v[np.maximum(idx, 0)] - self.v[0], 0.0)


class PDLVECM:
    """Microstructure VAR/VECM with PDL lags (Hasbrouck's MVARi), estimated from sparse data."""

    def __init__(self, prices: Sequence[PriceSeries], polys: Sequence[Poly], second_factor: int,
                 ecm: bool = True, intercept: bool = False):
        self.prices = list(prices)
        self.polys = list(polys)
        self.second_factor = int(second_factor)
        self.n = len(self.prices)
        self.ecm = ecm and self.n >= 2
        self.intercept = intercept
        self.max_lag = max(p.k_offset + p.n for p in self.polys)
        self.price_names = [p.name for p in self.prices]
        self._layout()
        self.first = self.last = None
        self.ecm_means = None
        self.zpz = None
        self.b = self.seb = self.tb = self.ecov = self.ecorr = None
        self.irft = self.irf = None

    # ------------------------------------------------------------------ layout
    def _layout(self):
        """Column names of the notional data matrix z (setNamesPDL):
        [const, (dP1, dP1p0d0, dP1p1d0, ...), (dP2, ...), ..., ecm terms]."""
        n_p = len(self.polys)
        self.block = 1 + n_p
        names = ["const"]
        for p in self.prices:
            names.append(f"d{p.name}")
            names += [f"d{p.name}{poly.var_name}" for poly in self.polys]
        self.ecm_names = [f"{self.prices[0].name}-{p.name}" for p in self.prices[1:]] if self.ecm else []
        names += self.ecm_names
        self.z_names = names
        self.lhs = [1 + j * self.block for j in range(self.n)]
        rhs = [0] if self.intercept else []
        for j in range(self.n):
            rhs += [1 + j * self.block + 1 + k for k in range(n_p)]
        rhs += [1 + self.n * self.block + m for m in range(len(self.ecm_names))]
        self.rhs = rhs
        self.y_names = [names[i] for i in self.lhs]
        self.x_names = [names[i] for i in self.rhs]

    # ------------------------------------------------------------- valid range
    def set_valid(self, first_time_s: float, last_time_s: float, log=print):
        """Estimation range in grid units (setValid). Hasbrouck's price objects have
        firstValid=1, so the lag constraint (firstValid + maxLag + 1) only binds at
        very coarse resolutions."""
        sf = self.second_factor
        self.first = max(int(round(first_time_s * sf)), 1 + self.max_lag + 1)
        self.last = int(round(last_time_s * sf))
        if self.last <= self.first:
            raise ValueError("empty estimation range")
        log(f"valid range: {self.first:12d} ({_hms(self.first / sf, sf)}) to "
            f"{self.last:12d} ({_hms(self.last / sf, sf)})")
        for p in self.prices:
            n_valid = int(np.sum((p.t >= self.first) & (p.t <= self.last)))
            log(f"{p.name:>12s} {len(p.t):10d} total observations, {n_valid:10d} in valid range.")
        if self.ecm:
            start = max(p.t[0] for p in self.prices)
            if start > self.first - 1:
                raise ValueError("a price series starts after the estimation range begins; "
                                 "error-correction terms are undefined")

    # ------------------------------------------------------------ cross products
    def _shifts(self) -> List[int]:
        s = {0, 1}
        for p in self.polys:
            s.update((p.k_offset + 1, p.k_offset + p.n + 1))
        return sorted(s)

    def _segments(self):
        """Breakpoints b_s in [first, last] and segment lengths w_s."""
        cands = [np.array([self.first], dtype=np.int64)]
        for p in self.prices:
            for s in self._shifts():
                c = p.t + s
                cands.append(c[(c > self.first) & (c <= self.last)])
        b = np.unique(np.concatenate(cands))
        w = np.diff(np.append(b, self.last + 1)).astype(np.float64)
        return b, w

    def _ecm_raw(self, q: np.ndarray) -> np.ndarray:
        """e_m(q-1) = p_1(q-1) - p_{m+1}(q-1), before demeaning."""
        p1 = self.prices[0].level(q - 1)
        out = np.column_stack([p1 - p.level(q - 1) for p in self.prices[1:]])
        if np.isnan(out).any():
            raise ValueError("error-correction term undefined inside the estimation range")
        return out

    def _z(self, q: np.ndarray) -> np.ndarray:
        z = np.empty((len(q), len(self.z_names)))
        z[:, 0] = 1.0
        for j, p in enumerate(self.prices):
            cache = {}

            def cum(shift):
                if shift not in cache:
                    cache[shift] = p.cum_dp(q - shift)
                return cache[shift]

            base = 1 + j * self.block
            z[:, base] = cum(0) - cum(1)
            for k, poly in enumerate(self.polys):
                z[:, base + 1 + k] = cum(poly.k_offset + 1) - cum(poly.k_offset + poly.n + 1)
        if self.ecm:
            z[:, 1 + self.n * self.block:] = self._ecm_raw(q) - self.ecm_means
        return z

    def build_zpz(self, chunk_rows: int = 500_000, log=print):
        """Exact cross-product matrix z'z over [first, last] (buildzpzSparsePDLpar).
        Error-correction terms are demeaned first (eVecDemean)."""
        if self.first is None:
            raise RuntimeError("call set_valid first")
        b, w = self._segments()
        n_obs = self.last - self.first + 1
        if self.ecm:
            tot = np.zeros(self.n - 1)
            for i in range(0, len(b), chunk_rows):
                tot += w[i:i + chunk_rows] @ self._ecm_raw(b[i:i + chunk_rows])
            self.ecm_means = tot / n_obs
            log("eVecMeans: " + "  ".join(f"{m:.5g}" for m in self.ecm_means))
        m = len(self.z_names)
        zpz = np.zeros((m, m))
        for i in range(0, len(b), chunk_rows):
            z = self._z(b[i:i + chunk_rows])
            zpz += z.T @ (z * w[i:i + chunk_rows, None])
        if zpz[0, 0] != n_obs:
            raise RuntimeError(f"segment weights sum to {zpz[0, 0]}, expected {n_obs}")
        self.zpz = zpz
        self.n_segments = len(b)
        return zpz

    # -------------------------------------------------------------- estimation
    def estimate(self):
        """OLS from the cross-product matrix (estimateSparse)."""
        n_obs = self.zpz[0, 0]
        xpx = self.zpz[np.ix_(self.rhs, self.rhs)]
        xpy = self.zpz[np.ix_(self.rhs, self.lhs)]
        ypy = self.zpz[np.ix_(self.lhs, self.lhs)]
        xpxi = np.linalg.inv(xpx)
        b = np.linalg.solve(xpx, xpy)
        ecov = (ypy - xpy.T @ b - b.T @ xpy + b.T @ xpx @ b) / n_obs
        ecov = (ecov + ecov.T) / 2
        v = np.diag(ecov)
        s = np.diag(1 / np.sqrt(v))
        self.b, self.ecov, self.ecorr = b, ecov, s @ ecov @ s
        self.seb = np.sqrt(np.outer(np.diag(xpxi), v))
        self.tb = b / self.seb
        return b

    def estimates_table(self) -> pd.DataFrame:
        """Coefficients and t statistics side by side, one column pair per equation."""
        cols = {}
        for i, y in enumerate(self.y_names):
            cols[y] = self.b[:, i]
            cols[f"t({y})"] = self.tb[:, i]
        return pd.DataFrame(cols, index=self.x_names)

    def _offset(self) -> int:
        return 1 if self.intercept else 0

    def _gamma(self) -> np.ndarray:
        """Error-correction loadings, n x (n-1)."""
        return self.b[-(self.n - 1):, :].T

    def phi_sum(self) -> np.ndarray:
        """Sum over all lags of the VAR coefficient matrices (sum(phiSet, 3)).
        With degree-0 PDLs, each coefficient covers n_k lags."""
        n_p, off = len(self.polys), self._offset()
        ps = np.zeros((self.n, self.n))
        for j in range(self.n):
            for k, poly in enumerate(self.polys):
                ps[:, j] += poly.n * self.b[off + j * n_p + k, :]
        return ps

    def theta_sum(self) -> np.ndarray:
        """Sum of the VMA coefficients from the Granger representation (thetaSum)."""
        if not self.ecm:
            raise RuntimeError("theta_sum requires the error-correction specification")
        gamma = self._gamma()
        bmat = np.hstack([np.ones((self.n - 1, 1)), -np.eye(self.n - 1)])
        gamma_null = null_space(gamma.T)
        b_null = null_space(bmat)
        mid = gamma_null.T @ (np.eye(self.n) - self.phi_sum()) @ b_null
        return b_null @ np.linalg.inv(mid) @ gamma_null.T

    # ---------------------------------------------------------- impulse responses
    def irf_packed(self, n_ahead: int, log=None):
        """Cumulative price impulse responses on a log-spaced horizon grid (irfPDLpacked).

        Returns (irft, irf) with irf[i, h, s] = response of price i at horizon irft[h]
        (in grid units) to a unit shock in price s at time 0. All shocks are run at
        once with the PDL state recursion (fPDLpacked), keeping a ring buffer of the
        last max_lag+1 price changes."""
        n, n_p, off = self.n, len(self.polys), self._offset()
        m = math.ceil(math.log10(n_ahead))
        k = np.unique(np.concatenate([(np.arange(11)[:, None] * 10 ** np.arange(m + 1)).ravel(), [n_ahead]]))
        irft = k[k <= n_ahead].astype(np.int64)

        coef = np.ascontiguousarray(self.b[off:off + n * n_p, :].T)     # n x (n*n_p)
        gamma = self._gamma() if self.ecm else None
        const = self.b[0][:, None] if self.intercept else None
        ko = np.array([p.k_offset for p in self.polys], dtype=np.int64)
        kn = ko + np.array([p.n for p in self.polys], dtype=np.int64)
        ring = self.max_lag + 1
        hist = np.zeros((ring, n, n))        # dp(t) for the last `ring` periods (time, price, shock)
        state = np.zeros((n, n_p, n))        # PDL window sums (price, poly, shock)
        p = np.zeros((n, n))                 # cumulative price (price, shock)
        out = np.empty((n, len(irft), n))
        eye = np.eye(n)
        nxt = 0
        for t in range(n_ahead + 1):
            dp = coef @ state.reshape(n * n_p, n)
            if gamma is not None:
                dp += gamma @ (p[0] - p[1:])
            if const is not None:
                dp += const
            if t == 0:
                dp += eye
            p += dp
            hist[t % ring] = dp
            state += (hist[(t - ko) % ring] - hist[(t - kn) % ring]).transpose(1, 0, 2)
            if t == irft[nxt]:
                out[:, nxt, :] = p
                nxt += 1
                if log and t >= 10 ** 6 and t % 10 ** 6 == 0:
                    log(f"  irf horizon {t:,} of {n_ahead:,}")
        self.irft, self.irf = irft, out
        return irft, out


def null_space(a: np.ndarray, rtol: float = 1e-10) -> np.ndarray:
    """Orthonormal basis for the null space of a (MATLAB null)."""
    a = np.atleast_2d(a)
    _, s, vt = np.linalg.svd(a)
    tol = rtol * (s[0] if len(s) else 1.0) * max(a.shape)
    rank = int(np.sum(s > tol))
    return vt[rank:].T.conj()


def _hms(seconds: float, second_factor: int) -> str:
    prec = max(0, int(round(math.log10(second_factor))))
    h = int(seconds // 3600)
    mnt = int((seconds % 3600) // 60)
    s = seconds % 60
    if prec == 0:
        return f"{h:02d}:{mnt:02d}:{int(s):02d}"
    return f"{h:02d}:{mnt:02d}:{s:0{prec + 3}.{prec}f}"
