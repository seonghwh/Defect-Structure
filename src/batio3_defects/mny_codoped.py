from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import brentq


# --- Physical / model constants (from your notebook) ---
K_B_EV_PER_K = 8.61733e-5  # Boltzmann constant in eV/K
B_SITE_DENSITY_CM3 = 1.55e22  # used in your notebook for (1 - A/B)*1.55E22


@dataclass(frozen=True)
class ReactionConstants:
    """Reaction constants at a single temperature."""
    KR: float
    Ki: float
    KS: float
    KMn43: float
    KMn32: float


def reaction_constants(TK: float) -> ReactionConstants:
    """
    Reproduce constants exactly as in MnY_codoped.ipynb (cell 2).
    """
    k = K_B_EV_PER_K

    KR = 1.06e71 * np.exp(-5.69 / (k * TK))
    Ki = 8.55e44 * np.exp(-2.91 / (k * TK))
    KS = (3.4e105) * np.exp(-2.795 / (k * TK))

    # Hagemann & Hennings (as used in notebook)
    KMn43 = 2 * 2.0e22 * np.exp(-(3.12 - 1.28) / (k * TK))
    KMn32 = 0.5 * 2.0e22 * np.exp(-(3.12 - 1.87) / (k * TK))

    return ReactionConstants(KR=KR, Ki=Ki, KS=KS, KMn43=KMn43, KMn32=KMn32)


def _mn_partition(n: float, Mn_total: float, KMn43: float, KMn32: float) -> Tuple[float, float, float, float]:
    """
    Matches your notebook definitions:
      Mn1 = Mn_total / (KMn43/n + 1 + n/KMn32)
      Mn2 = Mn1 * n / KMn32
      Mn0 = Mn_total - Mn1 - Mn2
    Returns: (Mn0, Mn1, Mn2, charge = Mn1 + 2*Mn2)
    """
    Mn1 = Mn_total / (KMn43 / n + 1.0 + n / KMn32)
    Mn2 = Mn1 * n / KMn32
    Mn0 = Mn_total - Mn1 - Mn2
    charge = Mn1 + 2.0 * Mn2
    return Mn0, Mn1, Mn2, charge


def _bracket_root(g, umin: float, umax: float, npts: int = 801) -> Tuple[float, float]:
    """
    Find [u_lo, u_hi] where g(u) changes sign. u is log10(n).
    """
    us = np.linspace(umin, umax, npts)
    vals = np.array([g(u) for u in us], dtype=float)

    # Skip NaNs/infs safely
    for i in range(len(us) - 1):
        v0, v1 = vals[i], vals[i + 1]
        if not np.isfinite(v0) or not np.isfinite(v1):
            continue
        if v0 == 0.0:
            return us[i], us[i]
        if np.sign(v0) * np.sign(v1) < 0:
            return us[i], us[i + 1]

    raise ValueError(
        "Could not bracket root in log10(n). "
        "Try widening umin/umax or inspect neutrality function."
    )


def solve_log10_n(neutrality_n, umin: float = -20.0, umax: float = 30.0) -> float:
    """
    Robust brentq solve in log-space:
      u = log10(n), solve neutrality(10^u) = 0
    """
    def g(u: float) -> float:
        return float(neutrality_n(10.0 ** u))

    u_lo, u_hi = _bracket_root(g, umin, umax)
    if u_lo == u_hi:
        return 10.0 ** u_lo

    u_root = brentq(g, u_lo, u_hi, maxiter=200)
    return 10.0 ** u_root


def solve_equilibrium_mny(
    pO2_grid: np.ndarray,
    TK: float,
    ratio_AB: float,
    Mn_total: float,
    Y_total: float,
) -> pd.DataFrame:
    """
    Equilibrium (high-T) Mn+Y co-doped case.
    Reproduces notebook cell 5 (ratio=1.0 Schottky) and cells 7/8 (ratio<1 fixed VBa).
    """
    rc = reaction_constants(TK)
    KR, Ki, KS, KMn43, KMn32 = rc.KR, rc.Ki, rc.KS, rc.KMn43, rc.KMn32

    rows = []
    for y in pO2_grid:
        y = float(y)

        def neutrality(n: float) -> float:
            p = Ki / n
            VO = KR / (n**2 * np.sqrt(y))
            Mn0, Mn1, Mn2, mn_charge = _mn_partition(n, Mn_total, KMn43, KMn32)

            if abs(ratio_AB - 1.0) < 1e-12:
                # Full Schottky: VBa = VTi, and KS couples to KR (matches your notebook)
                VTi = np.sqrt(KS / (KR**3)) * (n**3) * (y ** (3.0 / 4.0))
                VBa = VTi
                neg_ionic = 2.0 * VBa + 4.0 * VTi  # = 6*VTi
            else:
                # Ti-rich: VBa fixed by non-stoichiometry (your notebook assumption)
                VBa = (1.0 - ratio_AB) * B_SITE_DENSITY_CM3
                VTi = 0.0
                neg_ionic = 2.0 * VBa

            # n + (ionic negative charge) + (Mn negative charge) - (p + 2VO + Y) = 0
            return (n + neg_ionic + mn_charge) - (p + 2.0 * VO + Y_total)

        n = solve_log10_n(neutrality, umin=-20.0, umax=30.0)

        p = Ki / n
        VO = KR / (n**2 * np.sqrt(y))
        Mn0, Mn1, Mn2, _ = _mn_partition(n, Mn_total, KMn43, KMn32)

        if abs(ratio_AB - 1.0) < 1e-12:
            VTi = np.sqrt(KS / (KR**3)) * (n**3) * (y ** (3.0 / 4.0))
            VBa = VTi
        else:
            VBa = (1.0 - ratio_AB) * B_SITE_DENSITY_CM3
            VTi = 0.0

        rows.append(
            dict(
                pO2=y,
                log10_pO2=np.log10(y),
                n=n,
                p=p,
                VO2=VO,
                VBa2=VBa,
                VTi4=VTi,
                Mn0=Mn0,
                Mn1=Mn1,
                Mn2=Mn2,
                YBa=Y_total,
                ratio_AB=ratio_AB,
                TK=TK,
            )
        )

    return pd.DataFrame(rows)


def solve_quenched_mny(
    pO2_grid: np.ndarray,
    TQK: float,
    ratio_AB: float,
    Mn_total_quench: float,
    Y_total_quench: float,
    frozen_eq: pd.DataFrame,
) -> pd.DataFrame:
    """
    Quenched (low-T) case with *partial* re-equilibration:
      - VBa and VTi are frozen from high-T equilibrium
      - oxygen vacancies VO re-equilibrate at TQK with pO2
      - electrons/holes + Mn redox re-equilibrate at TQK
    """
    rcQ = reaction_constants(TQK)
    KRQ, KiQ, KMn43Q, KMn32Q = rcQ.KR, rcQ.Ki, rcQ.KMn43, rcQ.KMn32

    # Frozen ionic defects (except VO)
    VBa_f = frozen_eq["VBa2"].to_numpy(dtype=float)
    VTi_f = frozen_eq["VTi4"].to_numpy(dtype=float)

    rows = []
    for i, y in enumerate(pO2_grid):
        y = float(y)

        def neutrality(n: float) -> float:
            p = KiQ / n
            VO = KRQ / (n**2 * np.sqrt(y))  # <-- VO equilibrates at TQK
            _, _, _, mn_charge = _mn_partition(n, Mn_total_quench, KMn43Q, KMn32Q)

            if abs(ratio_AB - 1.0) < 1e-12:
                neg_ionic = 2.0 * VBa_f[i] + 4.0 * VTi_f[i]  # frozen
            else:
                neg_ionic = 2.0 * VBa_f[i]  # frozen

            return (n + neg_ionic + mn_charge) - (p + 2.0 * VO + Y_total_quench)

        n = solve_log10_n(neutrality, umin=-30.0, umax=35.0)

        p = KiQ / n
        VO = KRQ / (n**2 * np.sqrt(y))  # <-- final equilibrated VO
        Mn0, Mn1, Mn2, _ = _mn_partition(n, Mn_total_quench, KMn43Q, KMn32Q)

        rows.append(
            dict(
                pO2=y,
                log10_pO2=np.log10(y),
                n=n,
                p=p,
                VO2=VO,          # <-- NOT frozen anymore
                VBa2=VBa_f[i],   # frozen
                VTi4=VTi_f[i],   # frozen
                Mn0=Mn0,
                Mn1=Mn1,
                Mn2=Mn2,
                YBa=Y_total_quench,
                ratio_AB=ratio_AB,
                TQK=TQK,
            )
        )

    return pd.DataFrame(rows)


def idx_nearest_pO2(pO2_grid: np.ndarray, target_atm: float = 0.21) -> int:
    return int(np.argmin(np.abs(pO2_grid - target_atm)))