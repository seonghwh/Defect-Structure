from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "src"))

from batio3_defects.solver import B_SITE_DENSITY_CM3, KRModel, solve_equilibrium, solve_quenched


TK = 1150 + 273
TQK = 500 + 273
AIR_PO2 = 0.21
PO2_MIN = 1e-20
PO2_MAX = 1e5
DEFAULT_POINTS = 401

MN_TOTAL_CM3 = 7.7e19
Y_TOTAL_CM3 = 7.7e19
BACKGROUND_ACCEPTOR_CM3 = 100.0e-6 * B_SITE_DENSITY_CM3


@dataclass(frozen=True)
class SampleConfig:
    sample_id: str
    panel: str
    plot_title: str
    table_label: str
    ratio_ab: float
    mn_total: float = 0.0
    y_total: float = 0.0
    acc_cm3: float = 0.0
    acc_charge: int = 1
    kr_model: KRModel = "undoped_effective"


SAMPLES: tuple[SampleConfig, ...] = (
    SampleConfig(
        sample_id="undoped_ab0999",
        panel="a",
        plot_title=r"Undoped BaTiO$_3$ (A/B = 0.999)",
        table_label="Undoped BaTiO3, A/B = 0.999",
        ratio_ab=0.999,
        acc_cm3=BACKGROUND_ACCEPTOR_CM3,
        kr_model="undoped_effective",
    ),
    SampleConfig(
        sample_id="mn_ab0994",
        panel="b",
        plot_title=r"0.5 mol% Mn-BaTiO$_3$ (A/B = 0.994)",
        table_label="0.5 mol% Mn-BaTiO3, A/B = 0.994",
        ratio_ab=0.994,
        mn_total=MN_TOTAL_CM3,
        kr_model="mn_effective",
    ),
    SampleConfig(
        sample_id="mn_ab0999",
        panel="c",
        plot_title=r"0.5 mol% Mn-BaTiO$_3$ (A/B = 0.999)",
        table_label="0.5 mol% Mn-BaTiO3, A/B = 0.999",
        ratio_ab=0.999,
        mn_total=MN_TOTAL_CM3,
        kr_model="mn_effective",
    ),
    SampleConfig(
        sample_id="mny_ab0999",
        panel="d",
        plot_title=r"0.5 mol% Mn+Y-BaTiO$_3$ (A/B = 0.999)",
        table_label="0.5 mol% Mn+Y-BaTiO3, A/B = 0.999",
        ratio_ab=0.999,
        mn_total=MN_TOTAL_CM3,
        y_total=Y_TOTAL_CM3,
        kr_model="mn_effective",
    ),
)


SPECIES_ORDER = ("p", "n", "VO2", "VBa2", "VTi4", "Mn0", "Mn1", "Mn2", "YBa", "Acc")

SPECIES_STYLE = {
    "p": dict(label=r"$h^{\bullet}$", color="#D55E00", linestyle="-", linewidth=2.1),
    "n": dict(label=r"$e^{\prime}$", color="#0072B2", linestyle="-", linewidth=2.1),
    "VO2": dict(label=r"$V_{\mathrm{O}}^{\bullet\bullet}$", color="#009E73", linestyle="-", linewidth=2.1),
    "VBa2": dict(label=r"$V_{\mathrm{Ba}}^{\prime\prime}$", color="#CC79A7", linestyle="--", linewidth=1.9),
    "VTi4": dict(label=r"$V_{\mathrm{Ti}}^{\prime\prime\prime\prime}$", color="#8C6D31", linestyle="--", linewidth=1.9),
    "Mn0": dict(label=r"$\mathrm{Mn}_{\mathrm{Ti}}^{\times}$", color="#4D4D4D", linestyle="-", linewidth=1.9),
    "Mn1": dict(label=r"$\mathrm{Mn}_{\mathrm{Ti}}^{\prime}$", color="#E69F00", linestyle="-.", linewidth=1.9),
    "Mn2": dict(label=r"$\mathrm{Mn}_{\mathrm{Ti}}^{\prime\prime}$", color="#56B4E9", linestyle=":", linewidth=2.2),
    "YBa": dict(label=r"$\mathrm{Y}_{\mathrm{Ba}}^{\bullet}$", color="#7E57C2", linestyle="--", linewidth=1.9),
    "Acc": dict(label=r"$A^{\prime}$", color="#222222", linestyle=":", linewidth=1.8),
}

AIR_TABLE_COLUMNS = (
    ("sample", "Composition"),
    ("ratio_AB", "A/B"),
    ("n", "n"),
    ("p", "p"),
    ("VO2", "V_O**"),
    ("VBa2", "V_Ba''"),
    ("VTi4", "V_Ti''''"),
    ("Mn0", "Mn_Ti^x"),
    ("Mn1", "Mn_Ti'"),
    ("Mn2", "Mn_Ti''"),
    ("YBa", "Y_Ba*"),
    ("Acc", "A'"),
)


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 8,
            "axes.labelsize": 9,
            "axes.titlesize": 9,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def positive_log10(values: pd.Series) -> np.ndarray:
    arr = values.to_numpy(dtype=float)
    logged = np.full_like(arr, np.nan, dtype=float)
    mask = arr > 0.0
    logged[mask] = np.log10(arr[mask])
    return logged


def solve_sample(sample: SampleConfig, pO2_grid: np.ndarray) -> pd.DataFrame:
    eq = solve_equilibrium(
        pO2_grid,
        TK,
        sample.ratio_ab,
        Mn_total=sample.mn_total,
        Y_total=sample.y_total,
        acc_cm3=sample.acc_cm3,
        acc_charge=sample.acc_charge,
        kr_model=sample.kr_model,
    )
    quench = solve_quenched(
        pO2_grid,
        TQK,
        sample.ratio_ab,
        frozen_eq=eq,
        Mn_total_quench=sample.mn_total,
        Y_total_quench=sample.y_total,
        acc_cm3=sample.acc_cm3,
        acc_charge=sample.acc_charge,
        vo_equilibrates=True,
        kr_model=sample.kr_model,
    )
    quench.insert(0, "sample_id", sample.sample_id)
    quench.insert(1, "panel", sample.panel)
    quench.insert(2, "sample", sample.table_label)
    quench.insert(3, "sample_title", sample.plot_title)
    return quench


def build_figure_data(points: int) -> pd.DataFrame:
    pO2_grid = np.logspace(np.log10(PO2_MIN), np.log10(PO2_MAX), points)
    pO2_grid = np.unique(np.r_[pO2_grid, AIR_PO2])
    pO2_grid.sort()
    return pd.concat([solve_sample(sample, pO2_grid) for sample in SAMPLES], ignore_index=True)


def build_air_table() -> pd.DataFrame:
    pO2_air = np.array([AIR_PO2], dtype=float)
    rows: list[dict[str, object]] = []

    for sample in SAMPLES:
        quench = solve_sample(sample, pO2_air)
        row = quench.iloc[0]
        rows.append(
            {
                "sample_id": sample.sample_id,
                "sample": sample.table_label,
                "ratio_AB": sample.ratio_ab,
                "temperature_C": TQK - 273,
                "pO2_atm": AIR_PO2,
                "KR_model": row.KR_model,
                "n": row.n,
                "p": row.p,
                "VO2": row.VO2,
                "VBa2": row.VBa2,
                "VTi4": row.VTi4,
                "Mn0": row.Mn0,
                "Mn1": row.Mn1,
                "Mn2": row.Mn2,
                "YBa": row.YBa,
                "Acc": row.Acc,
            }
        )

    return pd.DataFrame(rows)


def species_present(df: pd.DataFrame, species: str) -> bool:
    if species not in df.columns:
        return False
    values = df[species].to_numpy(dtype=float)
    return bool(np.any(values > 0.0))


def add_air_marker(ax: plt.Axes, air_po2: float) -> None:
    air_x = np.log10(air_po2)
    ax.axvline(air_x, color="#B23A48", linestyle=(0, (4, 3)), linewidth=1.0, alpha=0.95)
    ax.text(
        air_x + 0.25,
        0.965,
        "air",
        transform=ax.get_xaxis_transform(),
        color="#B23A48",
        ha="left",
        va="top",
        fontsize=7.5,
    )


def plot_figure(
    data: pd.DataFrame,
    out_dir: Path,
    output_stem: str,
    formats: list[str],
    dpi: int,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> list[Path]:
    configure_matplotlib()
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.9), sharex=True, sharey=True)
    axes_flat = axes.ravel()
    used_species: set[str] = set()

    for ax, sample in zip(axes_flat, SAMPLES):
        df = data[data["sample_id"] == sample.sample_id]
        for species in SPECIES_ORDER:
            if not species_present(df, species):
                continue
            style = SPECIES_STYLE[species]
            ax.plot(
                df["log10_pO2"],
                positive_log10(df[species]),
                color=style["color"],
                linestyle=style["linestyle"],
                linewidth=style["linewidth"],
                solid_capstyle="round",
            )
            used_species.add(species)

        add_air_marker(ax, AIR_PO2)
        ax.set_title(f"({sample.panel}) {sample.plot_title}", loc="left", pad=5)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.xaxis.set_major_locator(MultipleLocator(5))
        ax.xaxis.set_minor_locator(MultipleLocator(1))
        ax.yaxis.set_major_locator(MultipleLocator(4))
        ax.yaxis.set_minor_locator(MultipleLocator(1))
        ax.grid(True, which="major", color="#D7DCE2", linewidth=0.6, alpha=0.9)
        ax.grid(True, which="minor", color="#ECEFF3", linewidth=0.4, alpha=0.7)
        ax.tick_params(direction="out", length=3.2, width=0.75)
        ax.tick_params(which="minor", direction="out", length=1.8, width=0.55)

    fig.supxlabel(r"$\log_{10}(p_{\mathrm{O}_2}/\mathrm{atm})$", y=0.108)
    fig.supylabel(r"$\log_{10}(c/\mathrm{cm}^{-3})$", x=0.035)

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=SPECIES_STYLE[species]["color"],
            linestyle=SPECIES_STYLE[species]["linestyle"],
            linewidth=SPECIES_STYLE[species]["linewidth"],
            label=SPECIES_STYLE[species]["label"],
        )
        for species in SPECIES_ORDER
        if species in used_species
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=5,
        frameon=False,
        bbox_to_anchor=(0.5, 0.015),
        handlelength=2.5,
        columnspacing=1.4,
        handletextpad=0.45,
    )
    fig.subplots_adjust(left=0.095, right=0.99, top=0.945, bottom=0.19, wspace=0.12, hspace=0.25)

    outputs: list[Path] = []
    for fmt in formats:
        fmt_clean = fmt.lower().lstrip(".")
        out_path = out_dir / f"{output_stem}.{fmt_clean}"
        fig.savefig(out_path, dpi=dpi, bbox_inches="tight", pad_inches=0.03)
        outputs.append(out_path)
    plt.close(fig)
    return outputs


def fmt_value(value: object) -> str:
    if isinstance(value, str):
        return value
    numeric = float(value)
    if numeric == 0.0:
        return "0"
    return f"{numeric:.3e}"


def write_air_table_markdown(air_table: pd.DataFrame, out_path: Path) -> None:
    headers = [label for _, label in AIR_TABLE_COLUMNS]
    align = [":--"] + ["--:"] * (len(headers) - 1)
    lines = [
        "# Table 3 Air Values",
        "",
        "Calculated at T = 500 C and pO2 = 0.21 atm after freezing cation vacancies from the 1150 C equilibrium state.",
        "Concentrations are in cm^-3.",
        "",
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(align) + " |",
    ]

    for _, row in air_table.iterrows():
        formatted: list[str] = []
        for key, _ in AIR_TABLE_COLUMNS:
            if key == "sample":
                formatted.append(str(row[key]))
            elif key == "ratio_AB":
                formatted.append(f"{float(row[key]):.3f}")
            else:
                formatted.append(fmt_value(row[key]))
        lines.append("| " + " | ".join(formatted) + " |")

    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate manuscript Figure 6 and the exact air-point defect concentration table."
    )
    parser.add_argument("--figures-dir", type=Path, default=REPO_ROOT / "figures")
    parser.add_argument("--results-dir", type=Path, default=REPO_ROOT / "results")
    parser.add_argument("--figure-stem", default="figure6_defect_concentrations")
    parser.add_argument("--table-stem", default="table3_air_values")
    parser.add_argument("--formats", nargs="+", default=["png", "pdf", "svg"])
    parser.add_argument("--dpi", type=int, default=600)
    parser.add_argument("--points", type=int, default=DEFAULT_POINTS)
    parser.add_argument("--xlim", nargs=2, type=float, default=[-20.0, 5.0], metavar=("MIN", "MAX"))
    parser.add_argument("--ylim", nargs=2, type=float, default=[5.0, 21.5], metavar=("MIN", "MAX"))
    parser.add_argument("--skip-figure", action="store_true")
    parser.add_argument("--skip-figure-data", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.results_dir.mkdir(parents=True, exist_ok=True)

    air_table = build_air_table()
    air_csv = args.results_dir / f"{args.table_stem}.csv"
    air_md = args.results_dir / f"{args.table_stem}.md"
    air_table.to_csv(air_csv, index=False)
    write_air_table_markdown(air_table, air_md)
    print(f"Saved air table: {air_csv}")
    print(f"Saved air table: {air_md}")

    if args.skip_figure:
        return

    figure_data = build_figure_data(points=args.points)
    if not args.skip_figure_data:
        data_out = args.results_dir / f"{args.figure_stem}_data.csv"
        figure_data.to_csv(data_out, index=False)
        print(f"Saved figure data: {data_out}")

    outputs = plot_figure(
        data=figure_data,
        out_dir=args.figures_dir,
        output_stem=args.figure_stem,
        formats=args.formats,
        dpi=args.dpi,
        xlim=(args.xlim[0], args.xlim[1]),
        ylim=(args.ylim[0], args.ylim[1]),
    )
    for out in outputs:
        print(f"Saved figure: {out}")


if __name__ == "__main__":
    main()
