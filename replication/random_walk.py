"""
Random-walk decomposition and information shares (port of randomWalkDecomp.m).

Given the VMA sum psi (one row of theta_sum) and the residual covariance Omega,
the random-walk variance is psi Omega psi'. For each ordering of the prices, the
Cholesky factor of Omega attributes that variance to the individual prices; the
information-share bounds are the min/max of each price's share across all n!
orderings. Grouped bounds sum the shares of a group's members within each
ordering before taking the min/max.
"""
from itertools import permutations
from typing import Dict, Sequence

import numpy as np
import pandas as pd

TRADING_SECONDS_PER_YEAR = 250 * 6.5 * 3600


class RandomWalkDecomp:
    def __init__(self, vma1: np.ndarray, ecov: np.ndarray, names: Sequence[str],
                 second_factor: int = None):
        self.vma1 = np.atleast_2d(vma1).mean(axis=0)   # rows should be identical
        self.ecov = np.asarray(ecov)
        self.names = list(names)
        self.n = len(self.names)
        self.second_factor = second_factor
        self.var_w = float(self.vma1 @ self.ecov @ self.vma1)
        self.var_w_annual = (self.var_w * TRADING_SECONDS_PER_YEAR * second_factor
                             if second_factor else None)
        self._all_perms()

    def _all_perms(self):
        perms = list(permutations(range(self.n)))
        vw = np.empty((len(perms), self.n))
        for r, k in enumerate(perms):
            k = list(k)
            chol = np.linalg.cholesky(self.ecov[np.ix_(k, k)])
            vw[r, k] = (self.vma1[k] @ chol) ** 2
        self.perms = perms
        self.vw_all_perms = vw
        self.is_all_perms = vw / vw.sum(axis=1).mean()
        self.is_bounds = pd.DataFrame({"min": self.is_all_perms.min(axis=0),
                                       "max": self.is_all_perms.max(axis=0)}, index=self.names)

    def grouped_bounds(self, groups: Dict[str, Sequence[str]]) -> pd.DataFrame:
        """Information-share bounds for groups of prices (isBoundsGrouped)."""
        rows = {}
        for label, members in groups.items():
            missing = [m for m in members if m not in self.names]
            if missing:
                raise KeyError(f"{missing} not among {self.names}")
            cols = [self.names.index(m) for m in members]
            share = self.is_all_perms[:, cols].sum(axis=1)
            rows[label] = (share.min(), share.max())
        self.is_group_bounds = pd.DataFrame.from_dict(rows, orient="index", columns=["min", "max"])
        return self.is_group_bounds

    def report(self) -> str:
        """Text block in the format of rDisplay."""
        lines = [f"Random walk decomposition for secondFactor={self.second_factor}",
                 "Sum of vma coefficients:",
                 pd.DataFrame([self.vma1], columns=self.names).to_string(index=False,
                                                                         float_format=lambda x: f"{x:.4g}"),
                 f"per period var_w:{self.var_w:12.6f} sd_w:{np.sqrt(self.var_w):12.6f}"]
        if self.var_w_annual is not None:
            lines.append(f"annualized var_w:{self.var_w_annual:12.6f} sd_w:{np.sqrt(self.var_w_annual):12.6f}")
        lines += ["info share bounds:", self.is_bounds.to_string(float_format=lambda x: f"{x:8.5f}")]
        if getattr(self, "is_group_bounds", None) is not None:
            lines += ["Grouped info share bounds",
                      self.is_group_bounds.to_string(float_format=lambda x: f"{x:8.5f}")]
        return "\n".join(lines)
