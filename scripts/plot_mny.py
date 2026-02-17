from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "results"
FIG_DIR = REPO_ROOT / "figures"


# Pretty labels (matches the style in your example)
LABELS = {
    "n": r"$n$",
    "p": r"$p$",
    "VO2": r"$V_{\mathrm{O}}^{\bullet\bullet}$",
    "VBa2": r"$V_{\mathrm{Ba}}^{\prime\prime}$",
    "VTi4": r"$V_{\mathrm{Ti}}^{''''}$",
    "Mn0": r"$\mathrm{Mn}_{\mathrm{Ti}}^{x}$",
    "Mn1": r"$\mathrm{Mn}_{\mathrm{Ti}}^{\prime}$",
    "Mn2": r"$\mathrm{Mn}_{\mathrm{Ti}}^{\prime\prime}$",
    "YBa": r"$\mathrm{Y}_{\mathrm{Ba}}^{\bullet}$",
}


def log10_safe(x: np.ndarray) -> np.ndarray:
    """log10 with zeros -> NaN (so they don't break the plot)."""
    x = np.asarray(x, dtype=float)
    out = np.full_like(x, np.nan, dtype=float)
    mask = x > 0
    out[mask] = np.log10(x[mask])
    return out


def plot_df(df: pd.DataFrame, title: str, outpath: Path, cols: list[str]) -> None:
    FIG_DIR.mkdir(exist_ok=True)

    x = df["log10_pO2"].to_numpy(dtype=float)
    air_x = np.log10(0.21)

    fig, ax = plt.subplots(figsize=(8, 5))

    for c in cols:
        if c not in df.columns:
            continue
        y = log10_safe(df[c].to_numpy(dtype=float))
        # Skip curves that are all NaN (e.g., VTi4=0 for A/B<1)
        if np.all(np.isnan(y)):
            continue
        ax.plot(x, y, label=LABELS.get(c, c))

    ax.axvline(air_x, linestyle="--", linewidth=1)

    ax.set_xlabel(r"$\log_{10}(p\mathrm{O}_2)\;(\mathrm{atm})$")
    ax.set_ylabel(r"$\log_{10}([\,]\;)\;(\mathrm{cm}^{-3})$")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", choices=["eq", "quench"], required=True, help="eq or quench CSV")
    parser.add_argument("--ratio", type=float, required=True, help="A/B ratio, e.g. 0.999")
    parser.add_argument("--title", type=str, default="", help="Optional custom title")
    args = parser.parse_args()

    state = args.state
    ratio = args.ratio

    csv = RESULTS_DIR / f"mny_{state}_ratio_{ratio:.3f}.csv"
    if not csv.exists():
        raise FileNotFoundError(f"Missing file: {csv}. Run scripts/run_mny_codoped.py first.")

    df = pd.read_csv(csv)

    # Use temperature column if present
    T = None
    if state == "eq" and "TK" in df.columns:
        T = int(df["TK"].iloc[0])
    if state == "quench" and "TQK" in df.columns:
        T = int(df["TQK"].iloc[0])

    # Default title similar to your example
    if args.title.strip():
        title = args.title.strip()
    else:
        # Adjust this text as you like
        Tc = f"{T - 273} °C" if T is not None else ""
        title = f"Mn+Y BTO (A/B = {ratio:.3f}); {state}; T = {Tc}".strip()

    # Columns to plot (edit if you want fewer/more)
    cols = ["n", "p", "VO2", "VBa2", "VTi4", "Mn0", "Mn1", "Mn2", "YBa"]

    out = FIG_DIR / f"mny_{state}_ratio_{ratio:.3f}.png"
    plot_df(df, title=title, outpath=out, cols=cols)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()