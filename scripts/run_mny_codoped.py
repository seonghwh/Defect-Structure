from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd

# Allow running without installing the package (simple dev workflow)
REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from batio3_defects.mny_codoped import (
    solve_equilibrium_mny,
    solve_quenched_mny,
    idx_nearest_pO2,
)

# ---- Defaults copied from your notebook ----
TK = 1150 + 273      # equilibrium temperature (K)
TQK = 500 + 273      # quenched temperature (K)
pO2_grid = np.logspace(-20, 5, 200)

ratios = [1.000, 0.999, 0.994]

# Equilibrium totals (your notebook used 0.5 mol% ~ 7.7E19 for Mn and Y)
Mn_eq = {r: 7.7e19 for r in ratios}
Y_eq = {r: 7.7e19 for r in ratios}

# Quenched totals (NOTE: your notebook is inconsistent: 1.000 & 0.999 use 2x, but 0.994 does not.)
# I’m preserving your notebook behavior exactly for now.
Mn_q = {r: 7.7e19 for r in ratios}
Y_q  = {r: 7.7e19 for r in ratios}


def main() -> None:
    results_dir = REPO_ROOT / "results"
    results_dir.mkdir(exist_ok=True)

    idx_air = idx_nearest_pO2(pO2_grid, 0.21)
    pO2_air = pO2_grid[idx_air]

    print(f"Using pO2_air ~ {pO2_air:.3g} atm (nearest grid point to 0.21 atm) at index {idx_air}\n")

    for r in ratios:
        # --- Equilibrium ---
        df_eq = solve_equilibrium_mny(
            pO2_grid=pO2_grid,
            TK=TK,
            ratio_AB=r,
            Mn_total=Mn_eq[r],
            Y_total=Y_eq[r],
        )
        out_eq = results_dir / f"mny_eq_ratio_{r:.3f}.csv"
        df_eq.to_csv(out_eq, index=False)

        # --- Quenched ---
        df_q = solve_quenched_mny(
            pO2_grid=pO2_grid,
            TQK=TQK,
            ratio_AB=r,
            Mn_total_quench=Mn_q[r],
            Y_total_quench=Y_q[r],
            frozen_eq=df_eq,
        )
        out_q = results_dir / f"mny_quench_ratio_{r:.3f}.csv"
        df_q.to_csv(out_q, index=False)

        # --- Print “air point” summary ---
        row_eq = df_eq.iloc[idx_air]
        row_q = df_q.iloc[idx_air]

        print(f"=== Mn+Y, A/B={r:.3f} ===")
        print(f"Equilibrium (TK={TK} K) @ pO2~{pO2_air:.3g} atm:")
        print(f"  n={row_eq.n:.2e}  p={row_eq.p:.2e}  VO2={row_eq.VO2:.2e}  VBa2={row_eq.VBa2:.2e}  VTi4={row_eq.VTi4:.2e}")
        print(f"  Mn0={row_eq.Mn0:.2e}  Mn1={row_eq.Mn1:.2e}  Mn2={row_eq.Mn2:.2e}  Y={row_eq.YBa:.2e}")

        print(f"Quenched (TQK={TQK} K) @ pO2~{pO2_air:.3g} atm (ionic frozen from eq):")
        print(f"  n={row_q.n:.2e}  p={row_q.p:.2e}  VO2={row_q.VO2:.2e}  VBa2={row_q.VBa2:.2e}  VTi4={row_q.VTi4:.2e}")
        print(f"  Mn0={row_q.Mn0:.2e}  Mn1={row_q.Mn1:.2e}  Mn2={row_q.Mn2:.2e}  Y={row_q.YBa:.2e}")
        print()

    print("Saved CSVs under: results/")


if __name__ == "__main__":
    main()