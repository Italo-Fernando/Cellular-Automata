from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

from analysis.analyze_results import summarize_daily_fractal_dimension, summarize_seir_dynamics
from box_counting import DEFAULT_BOX_SIZES, box_counting_curve
from cellular_automaton import STATE_I, VARIANT_PARAMETERS

# Paleta categórica de identidade das variantes (ordem fixa, validada para daltonismo).
VARIANT_COLORS = {"original": "#2a78d6", "delta": "#008300", "omicron": "#e87ba4"}
VARIANT_LINESTYLES = {"original": "-", "delta": "--", "omicron": ":"}

# Cores dos compartimentos SEIR, com mapeamento semântico
# (S = neutro, E = alerta, I = crítico, R = bom).
STATE_COLORS = {"S": "#898781", "E": "#fab219", "I": "#d03b3b", "R": "#0ca30c"}

TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRIDLINE = "#e1e0d9"
DPI = 300

plt.rcParams.update({
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": GRIDLINE,
    "axes.grid": True,
    "grid.color": GRIDLINE,
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "text.color": TEXT_PRIMARY,
    "axes.labelcolor": TEXT_PRIMARY,
    "xtick.color": TEXT_SECONDARY,
    "ytick.color": TEXT_SECONDARY,
    "font.size": 11,
})


def _variant_label(variant_name):
    return VARIANT_PARAMETERS[variant_name]["label"]


def _save(fig, output_path):
    fig.tight_layout()
    fig.savefig(output_path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def plot_grid_snapshots(representative_runs, snapshot_days=(5, 20, 50, 100), output_path=None):
    """Grade colorida por estado (S/E/I/R) em instantes selecionados, lado a lado por variante."""
    cmap = ListedColormap([STATE_COLORS["S"], STATE_COLORS["E"], STATE_COLORS["I"], STATE_COLORS["R"]])
    variants = list(representative_runs.keys())

    fig, axes = plt.subplots(len(variants), len(snapshot_days),
                              figsize=(2.3 * len(snapshot_days), 2.5 * len(variants)))
    axes = np.atleast_2d(axes)

    for row, variant_name in enumerate(variants):
        history = representative_runs[variant_name]
        for col, target_day in enumerate(snapshot_days):
            idx = min(target_day, len(history) - 1)
            actual_day, grid = history[idx]
            ax = axes[row, col]
            ax.imshow(grid, cmap=cmap, vmin=0, vmax=3, interpolation="nearest")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.grid(False)
            for spine in ax.spines.values():
                spine.set_visible(False)
            if row == 0:
                ax.set_title(f"t = {actual_day}", fontsize=10)
            if col == 0:
                ax.set_ylabel(_variant_label(variant_name), fontsize=10)

    legend_handles = [Patch(color=STATE_COLORS[state], label=state) for state in ("S", "E", "I", "R")]
    fig.legend(handles=legend_handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.03))
    _save(fig, output_path)
    return fig


def plot_epidemic_curves(results, output_path=None):
    seir_summary = summarize_seir_dynamics(results)
    variants = list(seir_summary["variant"].unique())

    fig, axes = plt.subplots(1, len(variants), figsize=(5 * len(variants), 4), sharey=True)
    axes = np.atleast_1d(axes)

    for ax, variant_name in zip(axes, variants):
        group = seir_summary[seir_summary["variant"] == variant_name]
        for col, state in (("n_s", "S"), ("n_e", "E"), ("n_i", "I"), ("n_r", "R")):
            ax.plot(group["day"], group[col], color=STATE_COLORS[state], linewidth=2, label=state)
        ax.set_title(_variant_label(variant_name))
        ax.set_xlabel("Tempo (dias)")

    axes[0].set_ylabel("Número de indivíduos")
    axes[0].legend(frameon=False)
    _save(fig, output_path)
    return fig


def plot_fractal_dimension_comparison(results, output_path=None):
    daily_summary = summarize_daily_fractal_dimension(results)

    fig, ax = plt.subplots(figsize=(7, 5))
    for variant_name, group in daily_summary.groupby("variant"):
        color = VARIANT_COLORS[variant_name]
        ax.plot(group["day"], group["mean"], color=color, linestyle=VARIANT_LINESTYLES[variant_name],
                linewidth=2, label=_variant_label(variant_name))
        ax.fill_between(group["day"], group["mean"] - group["std"], group["mean"] + group["std"],
                         color=color, alpha=0.15, linewidth=0)

    ax.set_xlabel("Tempo (dias)")
    ax.set_ylabel("Dimensão de box-counting D(t)")
    ax.legend(frameon=False, title="Variante")
    _save(fig, output_path)
    return fig


def plot_box_counting_example(binary_grid, box_sizes=DEFAULT_BOX_SIZES, output_path=None, title=None):
    sizes, counts = box_counting_curve(binary_grid, box_sizes)
    if len(sizes) < 2:
        return None

    log_inv_b = np.log(1.0 / sizes)
    log_n = np.log(counts)
    slope, intercept = np.polyfit(log_inv_b, log_n, 1)

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    ax.scatter(log_inv_b, log_n, color=VARIANT_COLORS["original"], s=50, zorder=3, label="N(b) observado")
    fit_x = np.linspace(log_inv_b.min(), log_inv_b.max(), 100)
    ax.plot(fit_x, slope * fit_x + intercept, color=TEXT_SECONDARY, linewidth=2,
            label=f"Ajuste linear (D = {slope:.3f})")

    ax.set_xlabel("log(1/b)")
    ax.set_ylabel("log N(b)")
    if title:
        ax.set_title(title)
    ax.legend(frameon=False)
    _save(fig, output_path)
    return fig


def generate_all_plots(results, representative_runs, output_dir):
    """Gera e salva os quatro gráficos exigidos pelo projeto em output_dir."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_grid_snapshots(representative_runs, output_path=output_dir / "grid_snapshots.png")
    plot_epidemic_curves(results, output_path=output_dir / "epidemic_curves.png")
    plot_fractal_dimension_comparison(results, output_path=output_dir / "fractal_dimension_comparison.png")

    # Ômicron forma clusters grandes e bem definidos cedo, o que rende uma
    # regressão log-log mais limpa para fins ilustrativos do que uma variante
    # com clusters pequenos/fragmentados.
    example_variant = "omicron" if "omicron" in representative_runs else next(iter(representative_runs))
    example_history = representative_runs[example_variant]
    example_day, example_grid = example_history[len(example_history) // 2]
    plot_box_counting_example(
        example_grid == STATE_I,
        output_path=output_dir / "box_counting_example.png",
        title=f"{_variant_label(example_variant)}, t = {example_day}",
    )
