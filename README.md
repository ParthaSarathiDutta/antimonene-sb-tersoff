# Antimonene Sb–Sb Tersoff Force Field

Training data, published LAMMPS Tersoff parameters (**ML-Tersoff**, **ML-Tersoff-1**), and a representative fitting workflow for antimonene.

## Contents

| Path | Description |
|------|-------------|
| [`parameters/ML-Tersoff_Sb.tersoff`](parameters/ML-Tersoff_Sb.tersoff) | ML-Tersoff |
| [`parameters/ML-Tersoff-1_Sb.tersoff`](parameters/ML-Tersoff-1_Sb.tersoff) | ML-Tersoff-1 |
| [`data/training/`](data/training/) | `1_beta.data`, `2_alpha.data` |
| [`fitting/`](fitting/) | Representative workflow code and `model.json` |

## LAMMPS

```lammps
pair_style tersoff
pair_coeff * * ML-Tersoff_Sb.tersoff Sb
```

```lammps
pair_coeff * * ML-Tersoff-1_Sb.tersoff Sb
```

## BLAST

The [BLAST](https://www.anl.gov/pse/blast-bridging-lengthtime-scales-using-atomistic-simulation-toolkit) framework developed at the [Center for Nanoscale Materials](https://www.anl.gov/cnm/about-the-cnm) at Argonne National Laboratory was used for fitting this model for antimonene. BLAST (Bridging Length/time scales via Atomistic Simulation Toolkit) is a multi-fidelity scale-bridging framework for training classical interatomic potentials for molecular simulations.

## Fitting workflow

The `fitting/` directory contains representative code describing the hierarchical property evaluation and optimization workflow used in this study.

- `evaluator.py` — representative hierarchical property evaluator
- `optimizer.py` — representative optimization interface; continuous MCTS was used in this study
- `model.json` — Tersoff parameter definitions and search ranges
- `settings.json` — example calculation settings

The evaluator and optimizer can be adapted to other fitting objectives and optimization algorithms.

## License

MIT ([`LICENSE`](LICENSE)) for scripts. CC BY 4.0 ([`LICENSE-CC-BY-4.0`](LICENSE-CC-BY-4.0)) for data and parameters.

## Citation

```bibtex
@software{antimonene_sb_tersoff,
  author  = {Dutta, Partha Sarathi and Koneru, Aditya and Muhammed, Adil and Manna, Sukriti and Chan, Henry and Loeffler, Troy and Sasikumar, Kiran and Sankaranarayanan, Subramanian KRS},
  title   = {Antimonene Sb--Sb Tersoff Force Field: Training Data, Parameters, and Representative Fitting Workflow},
  year    = {2026},
  url     = {https://github.com/ParthaSarathiDutta/antimonene-sb-tersoff},
  version = {v1.0.1}
}
```

Calculations used NERSC Perlmutter.
