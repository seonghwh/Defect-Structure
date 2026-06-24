from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

# Allow running without installing the package (simple dev workflow)
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from batio3_defects.solver import solve_equilibrium, solve_quenched

# ---- Defaults copied from your notebook ----
TK = 1150 + 273      # equilibrium temperature (K)
TQK = 500 + 273      # quenched temperature (K)
KR_MODEL = "mn_effective"
pO2_grid = np.logspace(-20, 5, 200)
pO2_air_exact = np.array([0.21])

ratios = [1.000, 0.999, 0.994]

# Equilibrium totals (your notebook used 0.5 mol% ~ 7.7E19 for Mn and Y)
Mn_eq = {r: 7.7e19 for r in ratios}
Y_eq = {r: 7.7e19 for r in ratios}

# Quenched totals
Mn_q = {r: 7.7e19 for r in ratios}
Y_q  = {r: 7.7e19 for r in ratios}


def main() -> None:
    results_dir = REPO_ROOT / "results"
    results_dir.mkdir(exist_ok=True)

    print(f"Using exact pO2_air = 0.21 atm for air-point summaries; KR_model={KR_MODEL}\n")

    for r in ratios:
        # Full-grid calculations for defect diagrams / plotting
        df_eq = solve_equilibrium(pO2_grid, TK, r, Mn_total=Mn_eq[r], Y_total=Y_eq[r], kr_model=KR_MODEL)
        df_q = solve_quenched(
            pO2_grid,
            TQK,
            r,
            frozen_eq=df_eq,
            Mn_total_quench=Mn_q[r],
            Y_total_quench=Y_q[r],
            vo_equilibrates=True,
            kr_model=KR_MODEL,
        )

        # Exact-air calculations for manuscript table / printed summaries
        df_eq_air = solve_equilibrium(pO2_air_exact, TK, r, Mn_total=Mn_eq[r], Y_total=Y_eq[r], kr_model=KR_MODEL)
        df_q_air = solve_quenched(
            pO2_air_exact,
            TQK,
            r,
            frozen_eq=df_eq_air,
            Mn_total_quench=Mn_q[r],
            Y_total_quench=Y_q[r],
            vo_equilibrates=True,
            kr_model=KR_MODEL,
        )

        out_q = results_dir / f"mny_quench_ratio_{r:.3f}.csv"
        df_q.to_csv(out_q, index=False)

        # Save exact-air values separately for easier Table 4 generation
        out_air = results_dir / f"mny_quench_air_ratio_{r:.3f}.csv"
        df_q_air.to_csv(out_air, index=False)

        row_eq = df_eq_air.iloc[0]
        row_q = df_q_air.iloc[0]

        print(f"=== Mn+Y, A/B={r:.3f} ===")
        print(f"Equilibrium (TK={TK} K) @ exact pO2=0.21 atm:")
        print(f"  n={row_eq.n:.2e}  p={row_eq.p:.2e}  VO2={row_eq.VO2:.2e}  VBa2={row_eq.VBa2:.2e}  VTi4={row_eq.VTi4:.2e}")
        print(f"  Mn0={row_eq.Mn0:.2e}  Mn1={row_eq.Mn1:.2e}  Mn2={row_eq.Mn2:.2e}  Y={row_eq.YBa:.2e}")

        print(f"Quenched (TQK={TQK} K) @ exact pO2=0.21 atm:")
        print(f"  n={row_q.n:.2e}  p={row_q.p:.2e}  VO2={row_q.VO2:.2e}  VBa2={row_q.VBa2:.2e}  VTi4={row_q.VTi4:.2e}")
        print(f"  Mn0={row_q.Mn0:.2e}  Mn1={row_q.Mn1:.2e}  Mn2={row_q.Mn2:.2e}  Y={row_q.YBa:.2e}")
        print()

    print("Saved CSVs under: results/")
    print("  - Full-grid quench outputs: mny_quench_ratio_*.csv")
    print("  - Exact-air quench outputs: mny_quench_air_ratio_*.csv")


if __name__ == "__main__":
    main()
