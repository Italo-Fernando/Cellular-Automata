# Cellula Automata

Modelagem matemático-computacional da propagação de variantes do SARS-CoV-2 com autômatos
celulares SEIR bidimensionais, medindo a dimensão fractal (box-counting de Kolmogorov) dos
clusters de infecção ao longo do tempo.

Projeto acadêmico: *Modelagem Matemático-Computacional da Propagação de Variantes do
SARS-CoV-2: Dimensão Fractal e Autômatos Celulares Aplicados ao Modelo SEIR*. Toda a simulação
é sintética (sem dados reais externos) — o objetivo é comparar como a trajetória D(t) da
dimensão fractal se comporta para três variantes (Original/Wuhan, Delta e Ômicron) com
perfis distintos de transmissibilidade, latência e período infeccioso.

## Como funciona

Uma grade fixa de 100×100 células (cada célula = um indivíduo) evolui em passos discretos de 1
dia, com vizinhança de Moore (8 vizinhos). Cada célula está em um dos quatro estados SEIR:

| Estado | Significado | Transição |
|---|---|---|
| **S** | Suscetível | vira **E** com prob. `1 − (1 − p)^k`, k = vizinhos em **I** |
| **E** | Exposto (infectado, ainda não transmite) | vira **I** com prob. `σ` |
| **I** | Infeccioso (transmite aos vizinhos) | vira **R** com prob. `γ` |
| **R** | Recuperado (imune permanente) | — |

A simulação começa com toda a grade em **S**, exceto a célula central em **I**, e termina quando
não resta nenhuma célula em **E** ou **I**. A cada dia, a dimensão fractal do cluster infeccioso
(estado **I**) é estimada por box-counting: cobre-se a grade com caixas de lados
`1, 2, 4, 8, 16, 32, 50`, conta-se `N(b)` (caixas com pelo menos uma célula infectada) e `D` é a
inclinação da reta de regressão de `log N(b)` contra `log(1/b)`.

### Parâmetros calibrados por variante

| Variante | p | σ (1/latência) | γ (1/infeccioso) | R0 |
|---|---|---|---|---|
| Original (Wuhan) | 0,04 | 0,15 | 0,11 | 2,79 |
| Delta | 0,06 | 0,23 | 0,10 | 5,08 |
| Ômicron | 0,23 | 0,29 | 0,19 | 9,5 |

## Estrutura do projeto

```
cellular_automaton.py     # classe SEIRCellularAutomaton + parâmetros das variantes
box_counting.py           # dimensão de box-counting (Kolmogorov)
simulation_runner.py      # ponto de entrada: roda as simulações e orquestra análise/gráficos
analysis/analyze_results.py  # métricas por execução, resumos e testes entre variantes
plots/generate_plots.py      # geração dos gráficos
output/                    # gerado ao rodar (CSVs e gráficos) — não versionado
```

## Instalação

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Como executar

```bash
python simulation_runner.py
```

Com os parâmetros padrão (30 execuções por variante, até 500 dias cada) a simulação completa
roda em cerca de 15-20 segundos em uma máquina comum.

### Flags disponíveis

| Flag | Padrão | O que faz |
|---|---|---|
| `--n-runs` | `30` | Número de repetições estocásticas (Monte Carlo) por variante. Mais execuções deixam a faixa de desvio padrão de D(t) mais estável, ao custo de mais tempo de execução. |
| `--seed` | `42` | Semente do gerador de números aleatórios, para reprodutibilidade. Sementes diferentes produzem trajetórias diferentes (a simulação é estocástica). |
| `--grid-size` | `100` | Lado da grade quadrada (número de células por lado). |
| `--max-days` | `500` | Número máximo de dias simulados por execução, mesmo que o surto ainda não tenha se extinguido (evita loops muito longos para variantes de espalhamento lento). |
| `--output-dir` | `output/` | Diretório onde os CSVs e a pasta `plots/` são salvos. |

Exemplo com uma rodada maior e uma pasta de saída diferente:

```bash
python simulation_runner.py --n-runs 50 --max-days 800 --output-dir output/rodada_final
```

## O que é gerado

Dentro de `--output-dir`:

- `simulation_results.csv` — série diária bruta (S, E, I, R e D(t)) de cada execução.
- `daily_fractal_dimension.csv` — média e desvio padrão de D(t) por variante e dia.
- `seir_dynamics_summary.csv` — contagens médias de S/E/I/R por variante e dia.
- `per_run_metrics.csv` — por execução: fração final do surto, pico de infectados (valor e dia),
  pico de D(t) (valor e dia) e inclinação de D(t) nos primeiros 20 dias.
- `summary_by_variant.csv` — as métricas acima resumidas (média ± desvio padrão) por variante.
- `pairwise_tests_peak_fractal_dim.csv` — teste de Mann-Whitney U entre cada par de variantes,
  sobre o pico de D(t).
- `plots/` — os quatro gráficos descritos abaixo, em PNG a 300 dpi.

## Gráficos gerados

Depois de rodar `simulation_runner.py`, os quatro gráficos ficam em `output/plots/`:

- **`grid_snapshots.png`** — estado da grade (S/E/I/R) em `t = 5, 20, 50, 100` dias, lado a lado
  por variante. Mostra visualmente por que a geometria dos clusters difere: a Ômicron forma uma
  frente de infecção compacta e quase circular, a Delta cresce de forma mais irregular, e a
  Original mal escapa da extinção estocástica no início do surto.
- **`epidemic_curves.png`** — contagem média de S, E, I, R ao longo do tempo, uma figura por
  variante (a curva epidêmica clássica).
- **`fractal_dimension_comparison.png`** — média ± desvio padrão de D(t) entre as execuções de
  Monte Carlo, as três variantes no mesmo gráfico. É o resultado central do trabalho: a Ômicron
  atinge D mais alto e mais cedo, mas cai abruptamente quando o cluster infeccioso se esgota
  (poucas células em **I** restantes); a Delta cresce mais devagar, mas também sobe mais e depois
  declina de forma semelhante quando seu próprio surto se aproxima da extinção; a Original
  permanece baixa e ruidosa, refletindo clusters pequenos e mais sujeitos à extinção estocástica.
- **`box_counting_example.png`** — ilustra, para um instante específico, como `D` é obtido a
  partir da reta de regressão de `log N(b)` contra `log(1/b)`.

## Resultados de exemplo

Tabela `summary_by_variant.csv` de uma rodada com 50 execuções por variante e até 800 dias
(`--n-runs 50 --max-days 800`):

| Variante | Fração final do surto | Pico de infectados (dia) | Pico de D(t) | Inclinação inicial de D(t) |
|---|---|---|---|---|
| Original (Wuhan) | 0,060 ± 0,063 | 21,1 (dia 252) | 0,53 ± 0,31 | 0,007 |
| Delta | 0,729 ± 0,414 | 245,3 (dia 327) | 0,85 ± 0,48 | 0,013 |
| Ômicron | 0,959 ± 0,198 | 420,9 (dia 179) | 1,16 ± 0,24 | 0,033 |

Testes de Mann-Whitney U par a par sobre o pico de D(t) confirmam que as três variantes diferem
entre si de forma estatisticamente significativa (p < 0,001 em todos os pares) nesta rodada.

Números variam de execução para execução (a simulação é estocástica); os arquivos CSV em
`output/` após rodar `simulation_runner.py` são a fonte de verdade para o texto do trabalho.

## Licença

Veja [LICENSE](LICENSE).
