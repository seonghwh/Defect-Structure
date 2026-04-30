from __future__ import annotations

from pathlib import Path
import sys
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from batio3_defects.solver import solve_equilibrium, solve_quenched, B_SITE_DENSITY_CM3

ACC_PPM = 100.0
ACC_CHARGE = 1  # A' (singly charged acceptor)
ACC_CM3 = ACC_PPM * 1e-6 * B_SITE_DENSITY_CM3   # ~1.55e18 cm^-3

TK = 1150 + 273
TQK = 500 + 273
pO2_grid = np.logspace(-20, 5, 200)
pO2_air_exact = np.array([0.21])
ratios = [1.000, 0.999, 0.994]


def main() -> None:
    results_dir = REPO_ROOT / "results"
    results_dir.mkdir(exist_ok=True)

    print("[undoped] Using exact pO2_air = 0.21 atm for air-point summaries")

    for r in ratios:
        # Full-grid calculations for defect diagrams / plotting
        df_eq = solve_equilibrium(
            pO2_grid,
            TK,
            r,
            Mn_total=0.0,
            Y_total=0.0,
            acc_cm3=ACC_CM3,
            acc_charge=ACC_CHARGE,
        )
        df_q = solve_quenched(
            pO2_grid,
            TQK,
            r,
            frozen_eq=df_eq,
            Mn_total_quench=0.0,
            Y_total_quench=0.0,
            acc_cm3=ACC_CM3,
            acc_charge=ACC_CHARGE,
            vo_equilibrates=True,
        )

        # Exact-air calculations for manuscript table / printed summaries
        df_eq_air = solve_equilibrium(
            pO2_air_exact,
            TK,
            r,
            Mn_total=0.0,
            Y_total=0.0,
            acc_cm3=ACC_CM3,
            acc_charge=ACC_CHARGE,
        )
        df_q_air = solve_quenched(
            pO2_air_exact,
            TQK,
            r,
            frozen_eq=df_eq_air,
            Mn_total_quench=0.0,
            Y_total_quench=0.0,
            acc_cm3=ACC_CM3,
            acc_charge=ACC_CHARGE,
            vo_equilibrates=True,
        )

        df_eq.to_csv(results_dir / f"undoped_eq_ratio_{r:.3f}.csv", index=False)
        df_q.to_csv(results_dir / f"undoped_quench_ratio_{r:.3f}.csv", index=False)
        df_q_air.to_csv(results_dir / f"undoped_quench_air_ratio_{r:.3f}.csv", index=False)

        row_eq = df_eq_air.iloc[0]
        row_q = df_q_air.iloc[0]
        print(f"  A/B={r:.3f} EQ @ exact pO2=0.21 atm: n={row_eq.n:.2e} p={row_eq.p:.2e} VO={row_eq.VO2:.2e}")
        print(f"          Q  @ exact pO2=0.21 atm: n={row_q.n:.2e} p={row_q.p:.2e} VO={row_q.VO2:.2e}")

    print("Saved CSVs under results/ (undoped_*)")
    print("  - Full-grid quench outputs: undoped_quench_ratio_*.csv")
    print("  - Exact-air quench outputs: undoped_quench_air_ratio_*.csv")
    print(f"[undoped] background acceptor: {ACC_PPM:g} ppm -> {ACC_CM3:.2e} cm^-3 (charge={ACC_CHARGE})")


if __name__ == "__main__":
    main()
