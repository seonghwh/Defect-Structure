from __future__ import annotations

import numpy as np
import pandas as pd

from batio3_defects.mny_codoped import reaction_constants, solve_log10_n, B_SITE_DENSITY_CM3


def _mn_partition(n: float, Mn_total: float, KMn43: float, KMn32: float):
    Mn1 = Mn_total / (KMn43 / n + 1.0 + n / KMn32)
    Mn2 = Mn1 * n / KMn32
    Mn0 = Mn_total - Mn1 - Mn2
    charge = Mn1 + 2.0 * Mn2
    return Mn0, Mn1, Mn2, charge


def solve_equilibrium_mn(
    pO2_grid: np.ndarray,
    TK: float,
    ratio_AB: float,
    Mn_total: float,
) -> pd.DataFrame:
    """
    Mn-doped BaTiO3 equilibrium at TK (no Y).
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
                VTi = np.sqrt(KS / (KR**3)) * (n**3) * (y ** (3.0 / 4.0))
                VBa = VTi
                neg_ionic = 2.0 * VBa + 4.0 * VTi
            else:
                VBa = (1.0 - ratio_AB) * B_SITE_DENSITY_CM3
                VTi = 0.0
                neg_ionic = 2.0 * VBa

            return (n + neg_ionic + mn_charge) - (p + 2.0 * VO)

        n = solve_log10_n(neutrality, umin=-30.0, umax=35.0)

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
                ratio_AB=ratio_AB,
                TK=TK,
                Mn_total=Mn_total,
            )
        )

    return pd.DataFrame(rows)


def solve_quenched_mn(
    pO2_grid: np.ndarray,
    TQK: float,
    ratio_AB: float,
    Mn_total_quench: float,
    frozen_eq: pd.DataFrame,
    vo_equilibrates: bool = True,
) -> pd.DataFrame:
    """
    Quenched at TQK:
      - VBa, VTi frozen from high-T eq
      - VO equilibrates at TQK by default (set vo_equilibrates=False to freeze it)
      - electrons/holes + Mn redox equilibrate at TQK
    """
    rcQ = reaction_constants(TQK)
    KRQ, KiQ, KMn43Q, KMn32Q = rcQ.KR, rcQ.Ki, rcQ.KMn43, rcQ.KMn32

    VBa_f = frozen_eq["VBa2"].to_numpy(dtype=float)
    VTi_f = frozen_eq["VTi4"].to_numpy(dtype=float)
    VO_f = frozen_eq["VO2"].to_numpy(dtype=float)

    rows = []
    for i, y in enumerate(pO2_grid):
        y = float(y)

        def neutrality(n: float) -> float:
            p = KiQ / n
            VO = (KRQ / (n**2 * np.sqrt(y))) if vo_equilibrates else VO_f[i]
            Mn0, Mn1, Mn2, mn_charge = _mn_partition(n, Mn_total_quench, KMn43Q, KMn32Q)

            if abs(ratio_AB - 1.0) < 1e-12:
                neg_ionic = 2.0 * VBa_f[i] + 4.0 * VTi_f[i]
            else:
                neg_ionic = 2.0 * VBa_f[i]

            return (n + neg_ionic + mn_charge) - (p + 2.0 * VO)

        n = solve_log10_n(neutrality, umin=-30.0, umax=35.0)

        p = KiQ / n
        VO = (KRQ / (n**2 * np.sqrt(y))) if vo_equilibrates else VO_f[i]
        Mn0, Mn1, Mn2, _ = _mn_partition(n, Mn_total_quench, KMn43Q, KMn32Q)

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
                ratio_AB=ratio_AB,
                TQK=TQK,
                Mn_total=Mn_total_quench,
            )
        )

    return pd.DataFrame(rows)