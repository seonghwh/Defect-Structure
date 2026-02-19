from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import brentq


# --- Shared constants (same as your current model) ---
K_B_EV_PER_K = 8.61733e-5
B_SITE_DENSITY_CM3 = 1.55e22  # used for (1 - A/B) * 1.55E22


@dataclass(frozen=True)
class ReactionConstants:
    KR: float
    Ki: float
    KS: float
    KMn43: float
    KMn32: float


def reaction_constants(TK: float) -> ReactionConstants:
    """Reproduce constants exactly as used in your notebook/model."""
    k = K_B_EV_PER_K
    KR = 1.06e71 * np.exp(-5.69 / (k * TK))
    Ki = 8.55e44 * np.exp(-2.91 / (k * TK))
    KS = 3.4e105 * np.exp(-2.795 / (k * TK))
    KMn43 = 2 * 2.0e22 * np.exp(-(3.12 - 1.28) / (k * TK))
    KMn32 = 0.5 * 2.0e22 * np.exp(-(3.12 - 1.87) / (k * TK))
    return ReactionConstants(KR=KR, Ki=Ki, KS=KS, KMn43=KMn43, KMn32=KMn32)


def idx_nearest_pO2(pO2_grid: np.ndarray, target_atm: float = 0.21) -> int:
    return int(np.argmin(np.abs(pO2_grid - target_atm)))


# --- Numerical solve utilities ---
def _bracket_root(g, umin: float, umax: float, npts: int = 1201) -> Tuple[float, float]:
    us = np.linspace(umin, umax, npts)
    vals = np.array([g(u) for u in us], dtype=float)
    for i in range(len(us) - 1):
        v0, v1 = vals[i], vals[i + 1]
        if not np.isfinite(v0) or not np.isfinite(v1):
            continue
        if v0 == 0.0:
            return us[i], us[i]
        if np.sign(v0) * np.sign(v1) < 0:
            return us[i], us[i + 1]
    raise ValueError("Could not bracket root in log10(n). Widen range or inspect neutrality.")


def solve_log10_n(neutrality_n, umin: float = -30.0, umax: float = 35.0) -> float:
    """Solve neutrality(n)=0 in log10(n) space for stability."""
    def g(u: float) -> float:
        return float(neutrality_n(10.0 ** u))

    u_lo, u_hi = _bracket_root(g, umin, umax)
    if u_lo == u_hi:
        return 10.0 ** u_lo
    u_root = brentq(g, u_lo, u_hi, maxiter=200)
    return 10.0 ** u_root


# --- Physics helpers ---
def mn_partition(n: float, Mn_total: float, KMn43: float, KMn32: float) -> Tuple[float, float, float, float]:
    """
    Mn1 = Mn_total / (KMn43/n + 1 + n/KMn32)
    Mn2 = Mn1 * n / KMn32
    Mn0 = Mn_total - Mn1 - Mn2
    charge = Mn1 + 2*Mn2
    """
    if Mn_total <= 0:
        return 0.0, 0.0, 0.0, 0.0
    Mn1 = Mn_total / (KMn43 / n + 1.0 + n / KMn32)
    Mn2 = Mn1 * n / KMn32
    Mn0 = Mn_total - Mn1 - Mn2
    charge = Mn1 + 2.0 * Mn2
    return Mn0, Mn1, Mn2, charge


def equilibrium_ionic_defects(
    n: float,
    pO2: float,
    ratio_AB: float,
    KR: float,
    KS: float,
) -> Tuple[float, float, float]:
    """
    Returns (VBa2, VTi4, neg_ionic_charge = 2*VBa + 4*VTi)

    - ratio_AB == 1.0: Schottky equilibrium (VBa = VTi from KS/KR relation)
    - ratio_AB < 1.0: Ti-rich (VBa fixed by non-stoichiometry, VTi=0)
    """
    if abs(ratio_AB - 1.0) < 1e-12:
        VTi = np.sqrt(KS / (KR**3)) * (n**3) * (pO2 ** (3.0 / 4.0))
        VBa = VTi
        neg_ionic = 2.0 * VBa + 4.0 * VTi
        return VBa, VTi, neg_ionic
    else:
        VBa = (1.0 - ratio_AB) * B_SITE_DENSITY_CM3
        VTi = 0.0
        neg_ionic = 2.0 * VBa
        return VBa, VTi, neg_ionic


# --- General solvers ---
def solve_equilibrium(
    pO2_grid: np.ndarray,
    TK: float,
    ratio_AB: float,
    Mn_total: float = 0.0,
    Y_total: float = 0.0,
    acc_cm3: float = 0.0,
    acc_charge: int = 1,
) -> pd.DataFrame:
    """
    General equilibrium solver:
      - Mn_total=0, Y_total=0 => undoped
      - Mn_total>0, Y_total=0 => Mn doped
      - Mn_total>0, Y_total>0 => Mn+Y codoped
    """
    rc = reaction_constants(TK)
    KR, Ki, KS, KMn43, KMn32 = rc.KR, rc.Ki, rc.KS, rc.KMn43, rc.KMn32
    acc_term = float(acc_charge) * float(acc_cm3)

    rows = []
    for y in pO2_grid:
        y = float(y)

        def neutrality(n: float) -> float:
            p = Ki / n
            VO = KR / (n**2 * np.sqrt(y))
            _, _, _, mn_charge = mn_partition(n, Mn_total, KMn43, KMn32)
            VBa, VTi, neg_ionic = equilibrium_ionic_defects(n, y, ratio_AB, KR, KS)
            return (n + neg_ionic + mn_charge + acc_term) - (p + 2.0 * VO + Y_total)

        n = solve_log10_n(neutrality)
        p = Ki / n
        VO = KR / (n**2 * np.sqrt(y))
        Mn0, Mn1, Mn2, _ = mn_partition(n, Mn_total, KMn43, KMn32)
        VBa, VTi, _ = equilibrium_ionic_defects(n, y, ratio_AB, KR, KS)

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
                Acc=acc_cm3,
                Acc_charge=acc_charge,
                ratio_AB=ratio_AB,
                TK=TK,
                Mn_total=Mn_total,
            )
        )

    return pd.DataFrame(rows)


def solve_quenched(
    pO2_grid: np.ndarray,
    TQK: float,
    ratio_AB: float,
    frozen_eq: pd.DataFrame,
    Mn_total_quench: float = 0.0,
    Y_total_quench: Optional[float] = None,
    acc_cm3: float = 0.0,
    acc_charge: int = 1,
    vo_equilibrates: bool = True,
) -> pd.DataFrame:
    """
    General quenched solver:
      - VBa and VTi are frozen from frozen_eq
      - VO is either equilibrated at TQK (vo_equilibrates=True) or frozen from frozen_eq
      - electrons/holes and Mn redox re-equilibrate at TQK
      - Y_total_quench defaults to frozen_eq["YBa"][0] if present, else 0.0
    """
    rcQ = reaction_constants(TQK)
    KRQ, KiQ, KMn43Q, KMn32Q = rcQ.KR, rcQ.Ki, rcQ.KMn43, rcQ.KMn32
    acc_term = float(acc_charge) * float(acc_cm3)

    VBa_f = frozen_eq["VBa2"].to_numpy(dtype=float)
    VTi_f = frozen_eq["VTi4"].to_numpy(dtype=float)
    VO_f = frozen_eq["VO2"].to_numpy(dtype=float)

    if Y_total_quench is None:
        if "YBa" in frozen_eq.columns:
            yvals = frozen_eq["YBa"].to_numpy(dtype=float)
            if not np.allclose(yvals, yvals[0]):
                raise ValueError("Frozen equilibrium YBa is not constant across pO2 grid.")
            Y_total_quench = float(yvals[0])
        else:
            Y_total_quench = 0.0

    rows = []
    for i, y in enumerate(pO2_grid):
        y = float(y)

        def neutrality(n: float) -> float:
            p = KiQ / n
            VO = (KRQ / (n**2 * np.sqrt(y))) if vo_equilibrates else VO_f[i]
            _, _, _, mn_charge = mn_partition(n, Mn_total_quench, KMn43Q, KMn32Q)

            if abs(ratio_AB - 1.0) < 1e-12:
                neg_ionic = 2.0 * VBa_f[i] + 4.0 * VTi_f[i]
            else:
                neg_ionic = 2.0 * VBa_f[i]

            return (n + neg_ionic + mn_charge + acc_term) - (p + 2.0 * VO + float(Y_total_quench))

        n = solve_log10_n(neutrality)
        p = KiQ / n
        VO = (KRQ / (n**2 * np.sqrt(y))) if vo_equilibrates else VO_f[i]
        Mn0, Mn1, Mn2, _ = mn_partition(n, Mn_total_quench, KMn43Q, KMn32Q)

        rows.append(
            dict(
                pO2=y,
                log10_pO2=np.log10(y),
                n=n,
                p=p,
                VO2=VO,
                VBa2=VBa_f[i],
                VTi4=VTi_f[i],
                Mn0=Mn0,
                Mn1=Mn1,
                Mn2=Mn2,
                YBa=float(Y_total_quench),
                Acc=acc_cm3,
                Acc_charge=acc_charge,
                ratio_AB=ratio_AB,
                TQK=TQK,
                Mn_total=Mn_total_quench,
                vo_equilibrates=bool(vo_equilibrates),
            )
        )

    return pd.DataFrame(rows)