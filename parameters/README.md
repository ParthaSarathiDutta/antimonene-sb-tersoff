# Published Tersoff parameters

Two Sb–Sb Tersoff parameter sets reported in the manuscript, fitted with the [BLAST](https://www.anl.gov/pse/blast-bridging-lengthtime-scales-using-atomistic-simulation-toolkit) framework ([CNM, Argonne National Laboratory](https://www.anl.gov/cnm/about-the-cnm)) against DFT reference properties for β- and α-antimonene.

## Files

| File | Manuscript name | Description |
|------|-----------------|-------------|
| [`ML-Tersoff_Sb.tersoff`](ML-Tersoff_Sb.tersoff) | **ML-Tersoff** | Best overall agreement with training and test properties |
| [`ML-Tersoff-1_Sb.tersoff`](ML-Tersoff-1_Sb.tersoff) | **ML-Tersoff-1** | Improved acoustic phonon branch scale; trade-off on elastic constants and test properties |

Both sets were selected from MCTS candidate parameter sets using primarily overall objective-function scores (lower is better), then evaluated for dynamical stability, thermal conductivity, mechanical fracture, and Grüneisen parameters. Reporting both illustrates trade-offs between structural, mechanical, vibrational, and thermal properties within the Tersoff functional form.

## ML-Tersoff (`ML-Tersoff_Sb.tersoff`)

Primary reported potential with the best overall agreement across fitted and test properties. Phonon dispersion shows slightly scaled-up acoustic branches relative to DFT.

Parameter line:

```
Sb Sb Sb 1 1.538065 3.383474 121310.341537 219.678679 -0.032348 2.037783 0.631679 0.708097 25.935178 3.763729 1.324866 3.281511 6077.106566
```

LAMMPS usage (metal units):

```lammps
pair_style tersoff
pair_coeff * * ML-Tersoff_Sb.tersoff Sb
```

## ML-Tersoff-1 (`ML-Tersoff-1_Sb.tersoff`)

Secondary reported potential with acoustic phonon branches closer in scale to the DFT reference, at the cost of poorer agreement for elastic constants and some test properties.

Parameter line:

```
Sb Sb Sb 1 1.728566 4.038711 473120.121738 475.785853 -0.064109 1.785647 0.887639 0.481402 19.620091 3.230913 1.187082 3.413096 12848.442002
```

LAMMPS usage (metal units):

```lammps
pair_style tersoff
pair_coeff * * ML-Tersoff-1_Sb.tersoff Sb
```

## Parameter format

`element element element m gamma lambda3 c d costheta0 n beta lambda2 B R D lambda1 A`

Ensure atom type 1 is assigned to Sb (mass 121.76 amu).

## License

Parameter files are released under CC BY 4.0.
