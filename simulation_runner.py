import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.analyze_results import run_full_analysis
from box_counting import DEFAULT_BOX_SIZES, box_counting_dimension
from cellular_automaton import (
    GRID_SIZE,
    MAX_DAYS,
    STATE_E,
    STATE_I,
    STATE_R,
    STATE_S,
    SEIRCellularAutomaton,
    VARIANT_PARAMETERS,
)
from plots.generate_plots import generate_all_plots

DEFAULT_N_RUNS = 30
OUTPUT_DIR = Path(__file__).resolve().parent / "output"


def _build_automaton(variant_name, grid_size, max_days, rng):
    params = VARIANT_PARAMETERS[variant_name]
    return SEIRCellularAutomaton(
        p=params["p"], sigma=params["sigma"], gamma=params["gamma"],
        grid_size=grid_size, max_days=max_days, rng=rng,
    )


def run_single_simulation(variant_name, run_id, seed=None, grid_size=GRID_SIZE,
                           max_days=MAX_DAYS, box_sizes=DEFAULT_BOX_SIZES):
    """Executa uma simulação e retorna a série diária de contagens SEIR e dimensão fractal."""
    rng = np.random.default_rng(seed)
    automaton = _build_automaton(variant_name, grid_size, max_days, rng)

    rows = []
    for day, grid in automaton.run():
        rows.append({
            "variant": variant_name,
            "run_id": run_id,
            "day": day,
            "n_s": int(np.sum(grid == STATE_S)),
            "n_e": int(np.sum(grid == STATE_E)),
            "n_i": int(np.sum(grid == STATE_I)),
            "n_r": int(np.sum(grid == STATE_R)),
            "fractal_dim_I": box_counting_dimension(grid == STATE_I, box_sizes),
        })
    return pd.DataFrame.from_records(rows)


def run_simulations(variants=None, n_runs=DEFAULT_N_RUNS, seed=42, grid_size=GRID_SIZE,
                     max_days=MAX_DAYS, box_sizes=DEFAULT_BOX_SIZES, output_dir=None):
    """Executa n_runs repetições de Monte Carlo por variante e concatena os resultados."""
    variants = variants or list(VARIANT_PARAMETERS.keys())
    master_rng = np.random.default_rng(seed)

    frames = []
    for variant_name in variants:
        for run_id in range(n_runs):
            run_seed = int(master_rng.integers(0, 2**32 - 1))
            frames.append(run_single_simulation(
                variant_name, run_id, seed=run_seed,
                grid_size=grid_size, max_days=max_days, box_sizes=box_sizes,
            ))

    results = pd.concat(frames, ignore_index=True)

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        results.to_csv(output_dir / "simulation_results.csv", index=False)

    return results


def run_representative_simulation(variant_name, seed=0, grid_size=GRID_SIZE, max_days=MAX_DAYS):
    """Executa uma única simulação guardando a grade completa a cada dia (para os snapshots espaciais)."""
    rng = np.random.default_rng(seed)
    automaton = _build_automaton(variant_name, grid_size, max_days, rng)
    return list(automaton.run())


def run_all_representative_simulations(variants=None, seed=0, grid_size=GRID_SIZE, max_days=MAX_DAYS):
    """Uma simulação representativa (com histórico completo de grades) por variante."""
    variants = variants or list(VARIANT_PARAMETERS.keys())
    return {
        variant_name: run_representative_simulation(
            variant_name, seed=seed + i, grid_size=grid_size, max_days=max_days,
        )
        for i, variant_name in enumerate(variants)
    }


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Simulação da propagação de variantes do SARS-CoV-2 via autômato celular SEIR "
                    "e análise da dimensão fractal dos clusters de infecção."
    )
    parser.add_argument("--n-runs", type=int, default=DEFAULT_N_RUNS,
                        help="Número de repetições de Monte Carlo por variante.")
    parser.add_argument("--seed", type=int, default=42, help="Semente para reprodutibilidade.")
    parser.add_argument("--grid-size", type=int, default=GRID_SIZE, help="Lado da grade quadrada.")
    parser.add_argument("--max-days", type=int, default=MAX_DAYS, help="Limite máximo de dias por simulação.")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR, help="Diretório de saída.")
    return parser.parse_args()


def main():
    args = _parse_args()
    variants = list(VARIANT_PARAMETERS.keys())

    print(f"Executando {args.n_runs} simulações por variante ({', '.join(variants)})...")
    results = run_simulations(
        variants=variants, n_runs=args.n_runs, seed=args.seed,
        grid_size=args.grid_size, max_days=args.max_days, output_dir=args.output_dir,
    )

    print("Executando simulações representativas para os snapshots espaciais...")
    representative_runs = run_all_representative_simulations(
        variants=variants, seed=args.seed, grid_size=args.grid_size, max_days=args.max_days,
    )

    print("Calculando estatísticas e testes entre variantes...")
    analysis = run_full_analysis(results, output_dir=args.output_dir)
    print(analysis["summary_by_variant"].to_string(index=False))
    print(analysis["pairwise_tests"].to_string(index=False))

    print("Gerando gráficos...")
    generate_all_plots(results, representative_runs, output_dir=args.output_dir / "plots")

    print(f"Concluído. Resultados e gráficos salvos em: {args.output_dir}")


if __name__ == "__main__":
    main()
