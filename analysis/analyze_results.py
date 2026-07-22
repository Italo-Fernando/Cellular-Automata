from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

EARLY_WINDOW_DAYS = 20


def summarize_daily_fractal_dimension(results, dim_col="fractal_dim_I"):
    """Média e desvio padrão de D(t) por variante e dia, entre execuções de Monte Carlo."""
    summary = results.groupby(["variant", "day"])[dim_col].agg(mean="mean", std="std").reset_index()
    summary["std"] = summary["std"].fillna(0.0)
    return summary


def summarize_seir_dynamics(results):
    return results.groupby(["variant", "day"])[["n_s", "n_e", "n_i", "n_r"]].mean().reset_index()


def summarize_per_run_metrics(results, early_window=EARLY_WINDOW_DAYS):
    """Métricas resumo por execução:

    - final_size_fraction: R(infinito) / população total
    - peak_infected / peak_infected_day: máximo de I(t) e quando ocorre
    - peak_fractal_dim / peak_fractal_day: máximo de D(t) e quando ocorre
    - early_slope_D: inclinação de D(t) nos primeiros `early_window` dias
    """
    rows = []
    for (variant, run_id), group in results.groupby(["variant", "run_id"]):
        group = group.sort_values("day")
        total_population = int(group.iloc[0][["n_s", "n_e", "n_i", "n_r"]].sum())
        final_size_fraction = float(group.iloc[-1]["n_r"] / total_population)

        peak_i_idx = group["n_i"].idxmax()
        peak_infected = int(group.loc[peak_i_idx, "n_i"])
        peak_infected_day = int(group.loc[peak_i_idx, "day"])

        valid_dim = group.dropna(subset=["fractal_dim_I"])
        if not valid_dim.empty:
            peak_d_idx = valid_dim["fractal_dim_I"].idxmax()
            peak_fractal_dim = float(valid_dim.loc[peak_d_idx, "fractal_dim_I"])
            peak_fractal_day = int(valid_dim.loc[peak_d_idx, "day"])
        else:
            peak_fractal_dim, peak_fractal_day = np.nan, np.nan

        early = valid_dim[valid_dim["day"] <= early_window]
        if len(early) >= 2:
            early_slope = float(np.polyfit(early["day"], early["fractal_dim_I"], 1)[0])
        else:
            early_slope = np.nan

        rows.append({
            "variant": variant,
            "run_id": run_id,
            "final_size_fraction": final_size_fraction,
            "peak_infected": peak_infected,
            "peak_infected_day": peak_infected_day,
            "peak_fractal_dim": peak_fractal_dim,
            "peak_fractal_day": peak_fractal_day,
            "early_slope_D": early_slope,
        })
    return pd.DataFrame(rows)


def summarize_by_variant(per_run_metrics):
    """Média e desvio padrão de cada métrica por variante, entre execuções — pronta para o texto."""
    return per_run_metrics.groupby("variant").agg(
        final_size_fraction_mean=("final_size_fraction", "mean"),
        final_size_fraction_std=("final_size_fraction", "std"),
        peak_infected_mean=("peak_infected", "mean"),
        peak_infected_day_mean=("peak_infected_day", "mean"),
        peak_fractal_dim_mean=("peak_fractal_dim", "mean"),
        peak_fractal_dim_std=("peak_fractal_dim", "std"),
        peak_fractal_day_mean=("peak_fractal_day", "mean"),
        early_slope_D_mean=("early_slope_D", "mean"),
        early_slope_D_std=("early_slope_D", "std"),
    ).reset_index()


def compare_variants_pairwise(per_run_metrics, metric_col):
    variants = sorted(per_run_metrics["variant"].unique())
    rows = []
    for variant_a, variant_b in combinations(variants, 2):
        sample_a = per_run_metrics.loc[per_run_metrics["variant"] == variant_a, metric_col].dropna()
        sample_b = per_run_metrics.loc[per_run_metrics["variant"] == variant_b, metric_col].dropna()

        if len(sample_a) > 0 and len(sample_b) > 0:
            u_stat, p_value = stats.mannwhitneyu(sample_a, sample_b, alternative="two-sided")
        else:
            u_stat, p_value = np.nan, np.nan

        rows.append({
            "metric": metric_col,
            "variant_a": variant_a,
            "variant_b": variant_b,
            "u_stat": u_stat,
            "p_value": p_value,
        })
    return pd.DataFrame(rows)


def run_full_analysis(results, output_dir=None):
    daily_summary = summarize_daily_fractal_dimension(results)
    seir_summary = summarize_seir_dynamics(results)
    per_run_metrics = summarize_per_run_metrics(results)
    summary_by_variant = summarize_by_variant(per_run_metrics)
    pairwise_tests = compare_variants_pairwise(per_run_metrics, "peak_fractal_dim")

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        daily_summary.to_csv(output_dir / "daily_fractal_dimension.csv", index=False)
        seir_summary.to_csv(output_dir / "seir_dynamics_summary.csv", index=False)
        per_run_metrics.to_csv(output_dir / "per_run_metrics.csv", index=False)
        summary_by_variant.to_csv(output_dir / "summary_by_variant.csv", index=False)
        pairwise_tests.to_csv(output_dir / "pairwise_tests_peak_fractal_dim.csv", index=False)

    return {
        "daily_summary": daily_summary,
        "seir_summary": seir_summary,
        "per_run_metrics": per_run_metrics,
        "summary_by_variant": summary_by_variant,
        "pairwise_tests": pairwise_tests,
    }
