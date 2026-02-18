from __future__ import annotations

import numpy as np
import pandas as pd

# Reuse tested utilities/constants from your MnY module
from batio3_defects.mny_codoped import reaction_constants, solve_log10_n, B_SITE_DENSITY_CM3


def solve_equilibrium_undoped(
    pO2_grid: np.ndarray,
    TK: float,
    ratio_AB: float,
    acc_cm3: float = 0.0,
    acc_charge: int = 1,
) -> pd.DataFrame:
    """
    Undoped BaTiO3 equilibrium at TK.
    Same canonical model as MnY but with Mn=0 and Y=0.
    """
    rc = reaction_constants(TK)
    KR, Ki, KS = rc.KR, rc.Ki, rc.KS

    rows = []
    for y in pO2_grid:
        y = float(y)

        def neutrality(n: float) -> float:
            p = Ki / n
            VO = KR / (n**2 * np.sqrt(y))

            if abs(ratio_AB - 1.0) < 1e-12:
                VTi = np.sqrt(KS / (KR**3)) * (n**3) * (y ** (3.0 / 4.0))
                VBa = VTi
                neg_ionic = 2.0 * VBa + 4.0 * VTi
            else:
                VBa = (1.0 - ratio_AB) * B_SITE_DENSITY_CM3
                VTi = 0.0
                neg_ionic = 2.0 * VBa

            acc = float(acc_charge) * float(acc_cm3)   # negative charge contribution
            return (n + neg_ionic + acc) - (p + 2.0 * VO)

        n = solve_log10_n(neutrality, umin=-30.0, umax=35.0)

        p = Ki / n
        VO = KR / (n**2 * np.sqrt(y))
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
                ratio_AB=ratio_AB,
                TK=TK,
                Acc=acc_cm3,
                Acc_charge=acc_charge,
            )
        )

    return pd.DataFrame(rows)


def solve_quenched_undoped(
    pO2_grid: np.ndarray,
    TQK: float,
    ratio_AB: float,
    frozen_eq: pd.DataFrame,
    vo_equilibrates: bool = True,
    acc_cm3: float = 0.0,
    acc_charge: int = 1,

) -> pd.DataFrame:
    """
    Quenched at TQK:
      - VBa, VTi frozen from high-T eq
      - VO either equilibrates at TQK (default) or is frozen (if vo_equilibrates=False)
      - electrons/holes equilibrate at TQK
    """
    rcQ = reaction_constants(TQK)
    KRQ, KiQ = rcQ.KR, rcQ.Ki

    VBa_f = frozen_eq["VBa2"].to_numpy(dtype=float)
    VTi_f = frozen_eq["VTi4"].to_numpy(dtype=float)
    VO_f = frozen_eq["VO2"].to_numpy(dtype=float)

    rows = []
    for i, y in enumerate(pO2_grid):
        y = float(y)

        def neutrality(n: float) -> float:
            p = KiQ / n
            VO = (KRQ / (n**2 * np.sqrt(y))) if vo_equilibrates else VO_f[i]

            if abs(ratio_AB - 1.0) < 1e-12:
                neg_ionic = 2.0 * VBa_f[i] + 4.0 * VTi_f[i]
            else:
                neg_ionic = 2.0 * VBa_f[i]

            acc = float(acc_charge) * float(acc_cm3)
            return (n + neg_ionic + acc) - (p + 2.0 * VO)

        n = solve_log10_n(neutrality, umin=-30.0, umax=35.0)
        p = KiQ / n
        VO = (KRQ / (n**2 * np.sqrt(y))) if vo_equilibrates else VO_f[i]

        rows.append(
            dict(
                pO2=y,
                log10_pO2=np.log10(y),
                n=n,
                p=p,
                VO2=VO,
                VBa2=VBa_f[i],
                VTi4=VTi_f[i],
                ratio_AB=ratio_AB,
                TQK=TQK,
                Acc=acc_cm3,
                Acc_charge=acc_charge,
            )
        )

    return pd.DataFrame(rows)