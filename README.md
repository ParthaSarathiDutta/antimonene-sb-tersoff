# Antimonene Sb–Sb Tersoff Force Field

Training data, published LAMMPS Tersoff parameters (**ML-Tersoff**, **ML-Tersoff-1**), and BLAST fitting scripts for antimonene.

## Contents

| Path | Description |
|------|-------------|
| [`parameters/ML-Tersoff_Sb.tersoff`](parameters/ML-Tersoff_Sb.tersoff) | ML-Tersoff |
| [`parameters/ML-Tersoff-1_Sb.tersoff`](parameters/ML-Tersoff-1_Sb.tersoff) | ML-Tersoff-1 |
| [`data/training/`](data/training/) | `1_beta.data`, `2_alpha.data` |
| [`fitting/`](fitting/) | BLAST scripts and `model.json` |

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

## Fitting (optional)

BLAST objective and search bounds are in [`fitting/main1.py`](fitting/main1.py) and [`fitting/model.json`](fitting/model.json). Requires [BLAST](https://www.anl.gov/pse/blast-bridging-lengthtime-scales-using-atomistic-simulation-toolkit) and LAMMPS.

## License

MIT (scripts). CC BY 4.0 (data and parameters).

## Citation

```bibtex
@software{antimonene_sb_tersoff,
  author = {Dutta, Partha Sarathi and Koneru, Aditya et al},
  title = {Antimonene Sb--Sb Tersoff Force Field: Training Data and BLAST Fitting Scripts},
  year = {2026},
  url = {https://github.com/ParthaSarathiDutta/antimonene-sb-tersoff}
}
```

Calculations used NERSC Perlmutter.
