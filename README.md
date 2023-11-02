# Dimension and degeneracy of zeros of parametric polynomial systems arising from reaction networks
This repository contains files for the manuscript [Dimension and degeneracy of zeros of parametric polynomial systems arising from reaction networks
](https://arxiv.org/abs/2304.02302) (2304.02302) by Elisenda Feliu, Oskar Henriksson, and Beatriz Pascual-Escudero.

## File descriptions
The repository contains the following files:
* A Maple file `functions.mpl` that contains functions for testing whether a network admits positive nondegenerate steady states when modeled with mass action kinetics.
* A file `investigated_models.csv` with all networks in ODEbase (as of November 2, 2023) with at least one reaction.
* A file `nondegenerate_networks` with all networks from `investigated_models.csv` that admit a positive nondegenerate steady state.
* A file `degenerate_networks` with all networks from `investigated_models.csv` that have a positive steady states, but all of them are degenerate.

## Example
Suppose we want to investigate the properties of a network with the following stoichiometric matrix and reactant matrix (this corresponds to the network `BIOMD0000000520` in ODEbase):

```
Gamma := Matrix([[-1, 0, 1, 0, 0, 0, 0], [0, 1, 0, -1, 0, 1, 0], [0, 0, 0, 0, 1, 0, -1]]);
B := Matrix([[1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 1, 1, 0], [0, 0, 0, 0, 0, 0, 1]]);
```

We begin by loading our Maple functions:

```
read("functions.mpl"):
```

The following command returns `true`, which means that the network admits positive steady states:

```
IsConsistent(Gamma)
```

The following command returns `false`, which means that all steady states are degenerate:

```
HasNondegenerateSteadyState(Gamma,B);
```

The following command returns `false`; with the notation from (4.2) in the paper, this means that all zeros of $f_\kappa$ are degenerate for all $\kappa$:
```
HasNondegenerateZeros(Gamma,B);
```
