# BaTiO3 Canonical Defect Model

This repository contains the reference Python implementation used to generate the canonical defect model results for the BaTiO3 manuscript on Mn doping, Mn+Y co-doping, and A/B cation ratio.

The current code output is the authoritative source for the manuscript defect-concentration values. Generated outputs are not tracked by git; reproduce them with the command below.

## Reproduce Figure 6 and Table 3

Create an environment and install the package:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
```

Generate the manuscript outputs:

```bash
python scripts/generate_manuscript_outputs.py
```

This writes:

- `figures/figure6_defect_concentrations.png`
- `figures/figure6_defect_concentrations.pdf`
- `figures/figure6_defect_concentrations.svg`
- `results/figure6_defect_concentrations_data.csv`
- `results/table3_air_values.csv`
- `results/table3_air_values.md`

Run the tests:

```bash
pytest
```

## Model Conditions

- High-temperature equilibration: 1150 C
- Low-temperature re-equilibration: 500 C
- Air point for Table 3: pO2 = 0.21 atm
- Site density used for A/B nonstoichiometry: `N_B = 1.55e22 cm^-3`
- Cation vacancies are frozen from the 1150 C equilibrium state.
- Electrons, holes, oxygen vacancies, and Mn charge states are re-equilibrated at 500 C.

## Current Table 3 Values

Calculated at 500 C and pO2 = 0.21 atm. Concentrations are in `cm^-3`.

| Composition | A/B | n | p | V_O** | V_Ba'' | V_Ti'''' | Mn_Ti^x | Mn_Ti' | Mn_Ti'' | Y_Ba* | A' |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| Undoped BaTiO3 | 0.999 | 4.188e+06 | 2.175e+19 | 5.401e+18 | 1.550e+19 | 0 | 0 | 0 | 0 | 0 | 1.550e+18 |
| 0.5 mol% Mn-BaTiO3 | 0.994 | 1.434e+07 | 6.351e+18 | 8.984e+19 | 9.300e+19 | 0 | 7.697e+19 | 3.421e+16 | 8.657e+09 | 0 | 0 |
| 0.5 mol% Mn-BaTiO3 | 0.999 | 3.598e+07 | 2.532e+18 | 1.428e+19 | 1.550e+19 | 0 | 7.691e+19 | 8.576e+16 | 5.444e+10 | 0 | 0 |
| 0.5 mol% Mn+Y-BaTiO3 | 0.999 | 4.775e+10 | 1.908e+15 | 8.106e+12 | 1.550e+19 | 0 | 3.104e+19 | 4.592e+19 | 3.869e+16 | 7.700e+19 | 0 |

## Repository Layout

- `src/batio3_defects/solver.py`: canonical defect model solver.
- `scripts/generate_manuscript_outputs.py`: only script needed to regenerate Figure 6 and Table 3 values.
- `tests/`: regression tests for the current air-point outputs and dopant site balances.
- `notebooks/archive/`: exploratory notebooks retained for provenance, not required for reproduction.
- `supporting_information/`: manuscript supporting-information drafts and model notes.
