from __future__ import annotations

from pathlib import Path
import sys
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from batio3_defects.undoped import solve_equilibrium_undoped, solve_quenched_undoped  # noqa: E402
from batio3_defects.mny_codoped import idx_nearest_pO2  # noqa: E402
from batio3_defects.mny_codoped import B_SITE_DENSITY_CM3

ACC_PPM = 100.0
ACC_CHARGE = 1  # A' (singly charged acceptor)
ACC_CM3 = ACC_PPM * 1e-6 * B_SITE_DENSITY_CM3   # ~1.55e18 cm^-3

TK = 1150 + 273
TQK = 500 + 273
pO2_grid = np.logspace(-20, 5, 200)
ratios = [1.000, 0.999, 0.994]


def main() -> None:
    results_dir = REPO_ROOT / "results"
    results_dir.mkdir(exist_ok=True)

    idx_air = idx_nearest_pO2(pO2_grid, 0.21)
    pO2_air = pO2_grid[idx_air]
    print(f"[undoped] Using pO2_air ~ {pO2_air:.3g} atm (index={idx_air})")

    for r in ratios:
        df_eq = solve_equilibrium_undoped(pO2_grid=pO2_grid, TK=TK, ratio_AB=r,
                                  acc_cm3=ACC_CM3, acc_charge=ACC_CHARGE)

        df_q = solve_quenched_undoped(pO2_grid=pO2_grid, TQK=TQK, ratio_AB=r, frozen_eq=df_eq,
                              vo_equilibrates=True, acc_cm3=ACC_CM3, acc_charge=ACC_CHARGE)
        
        df_eq.to_csv(results_dir / f"undoped_eq_ratio_{r:.3f}.csv", index=False)
        df_q.to_csv(results_dir / f"undoped_quench_ratio_{r:.3f}.csv", index=False)

        row_eq = df_eq.iloc[idx_air]
        row_q = df_q.iloc[idx_air]
        print(f"  A/B={r:.3f} EQ: n={row_eq.n:.2e} p={row_eq.p:.2e} VO={row_eq.VO2:.2e}")
        print(f"          Q : n={row_q.n:.2e} p={row_q.p:.2e} VO={row_q.VO2:.2e}")

    print("Saved CSVs under results/ (undoped_*)")
    print(f"[undoped] background acceptor: {ACC_PPM:g} ppm -> {ACC_CM3:.2e} cm^-3 (charge={ACC_CHARGE})")


if __name__ == "__main__":
    main()