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
Suppose we want to investigate the properties of a network with the following stoichiometric matrix and reactant matrix (this corresponds to `BIOMD0000000520` in ODEbase):

```
Gamma := Matrix([[-1, 0, 1, 0, 0, 0, 0], [0, 1, 0, -1, 0, 1, 0], [0, 0, 0, 0, 1, 0, -1]]);
B := Matrix([[1, 1, 1, 0, 0, 0, 0], [0, 0, 0, 1, 1, 1, 0], [0, 0, 0, 0, 0, 0, 1]]);
```

We begin by loading our Maple functions:

```
read("functions.mpl"):
```

In order to check whether there are positive steady states, we write:

```
IsConsistent(Gamma)
```
which, in this case, gives the output `true`.

In order to check whether there are nondegenerate steady states, we write:

```
HasNondegenerateSteadyState(Gamma,B);
```
which, in this case, gives the output `false`.

Finally, we check whether there are nondegenerate zeros of the steady state polynomials, by writing:
```
HasNondegenerateZeros(Gamma,B);
```
which, in this case, also gives the output `false`.
