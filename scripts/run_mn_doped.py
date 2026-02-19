from __future__ import annotations

from pathlib import Path
import sys
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from batio3_defects.solver import solve_equilibrium, solve_quenched, idx_nearest_pO2

TK = 1150 + 273
TQK = 500 + 273
pO2_grid = np.logspace(-20, 5, 200)
ratios = [1.000, 0.999, 0.994]

# Set your Mn dopant concentration here (cm^-3)
Mn_total = 7.7e19  # 0.5 mol% (same convention you used)


def main() -> None:
    results_dir = REPO_ROOT / "results"
    results_dir.mkdir(exist_ok=True)

    idx_air = idx_nearest_pO2(pO2_grid, 0.21)
    pO2_air = pO2_grid[idx_air]
    print(f"[mn] Using pO2_air ~ {pO2_air:.3g} atm (index={idx_air}); Mn_total={Mn_total:.3e} cm^-3")

    for r in ratios:
        df_eq = solve_equilibrium(pO2_grid, TK, r, Mn_total=Mn_total, Y_total=0.0)
        df_q  = solve_quenched(pO2_grid, TQK, r, frozen_eq=df_eq, Mn_total_quench=Mn_total, Y_total_quench=0.0, vo_equilibrates=True)

        df_eq.to_csv(results_dir / f"mn_eq_ratio_{r:.3f}.csv", index=False)
        df_q.to_csv(results_dir / f"mn_quench_ratio_{r:.3f}.csv", index=False)

        row_eq = df_eq.iloc[idx_air]
        row_q = df_q.iloc[idx_air]
        print(f"  A/B={r:.3f} EQ: n={row_eq.n:.2e} p={row_eq.p:.2e} VO={row_eq.VO2:.2e}")
        print(f"          Q : n={row_q.n:.2e} p={row_q.p:.2e} VO={row_q.VO2:.2e}")

    print("Saved CSVs under results/ (mn_*)")


if __name__ == "__main__":
    main()