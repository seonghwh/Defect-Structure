from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


REPO_ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = REPO_ROOT / "results"
FIG_DIR = REPO_ROOT / "figures"

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

def format_acceptor_note(df: pd.DataFrame) -> str:
    """
    Returns a short note like 'Acc=1.55e18 cm^-3 (charge=1)' if present and nonzero,
    otherwise returns ''.
    """
    if "Acc" not in df.columns:
        return ""
    acc = float(df["Acc"].iloc[0])
    if acc <= 0:
        return ""
    ch = int(df["Acc_charge"].iloc[0]) if "Acc_charge" in df.columns else 1
    return f"Acc={acc:.2e} cm^-3 (charge={ch})"

def log10_safe(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    out = np.full_like(x, np.nan, dtype=float)
    m = x > 0
    out[m] = np.log10(x[m])
    return out


def plot_df(df: pd.DataFrame, title: str, outpath: Path, cols: list[str], air_pO2: float = 0.21) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    x = df["log10_pO2"].to_numpy(dtype=float)
    air_x = np.log10(air_pO2)

    fig, ax = plt.subplots(figsize=(8, 5))

    for c in cols:
        if c not in df.columns:
            continue
        y = log10_safe(df[c].to_numpy(dtype=float))
        if np.all(np.isnan(y)):
            continue
        ax.plot(x, y, label=LABELS.get(c, c))

    # Air reference line (no label / no annotation)
    ax.axvline(np.log10(air_pO2), linestyle="--", linewidth=1)

    ax.set_xlabel(r"$\log_{10}(p\mathrm{O}_2)\;(\mathrm{atm})$")
    ax.set_ylabel(r"$\log_{10}([\,]\;)\;(\mathrm{cm}^{-3})$")
    ax.set_title(title)

    # Optional: show acceptor note on-plot for undoped
    if "Acc" in df.columns:
        acc = float(df["Acc"].iloc[0])
        if acc > 0:
            ch = int(df["Acc_charge"].iloc[0]) if "Acc_charge" in df.columns else 1
            ax.text(
                0.02, 0.98,
                f"Acc={acc:.2e} cm$^{{-3}}$ (charge={ch})",
                transform=ax.transAxes,
                ha="left", va="top",
            )
            
    ax.grid(True, alpha=0.3)

    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=True)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=["undoped", "mn", "mny"], required=True)
    ap.add_argument("--state", choices=["eq", "quench"], required=True)
    ap.add_argument("--ratio", type=float, required=True)
    ap.add_argument("--title", type=str, default="")
    ap.add_argument("--air-po2", type=float, default=0.21)
    args = ap.parse_args()

    csv = RESULTS_DIR / f"{args.case}_{args.state}_ratio_{args.ratio:.3f}.csv"
    if not csv.exists():
        raise FileNotFoundError(f"Missing file: {csv}. Run the corresponding run script first.")

    df = pd.read_csv(csv)
    acc_note = format_acceptor_note(df) if args.case == "undoped" else ""

    T = None
    if args.state == "eq" and "TK" in df.columns:
        T = int(df["TK"].iloc[0])
    if args.state == "quench" and "TQK" in df.columns:
        T = int(df["TQK"].iloc[0])

    if args.title.strip():
        title = args.title.strip()
    else:
        Tc = f"{T - 273} °C" if T is not None else ""
        title = f"{args.case.upper()} BTO (A/B={args.ratio:.3f}); {args.state}; T={Tc}".strip()
        if acc_note:
            title += f"; {acc_note}"

    cols = ["n", "p", "VO2", "VBa2", "VTi4", "Mn0", "Mn1", "Mn2", "YBa"]

    out = FIG_DIR / f"{args.case}_{args.state}_ratio_{args.ratio:.3f}.png"
    plot_df(df, title=title, outpath=out, cols=cols, air_pO2=args.air_po2)
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()