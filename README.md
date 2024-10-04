# Dimension and degeneracy of zeros of parametric polynomial systems arising from reaction networks
This repository contains files for the manuscript [Generic consistency and nondegeneracy of vertically parametrized systems
](https://arxiv.org/abs/2304.02302) (2304.02302) by Elisenda Feliu, Oskar Henriksson, and Beatriz Pascual-Escudero, as well 
as a forthcoming manuscript on applications to chemical reaction network theory.

> [!WARNING]  
> Parts of the code is still experimental. Use at own risk!

## File descriptions
The repository contains the following files:
* A Julia file `julia/functions.jl` that contains functions for testing whether a network admits positive nondegenerate steady states when modeled with generalized) mass action kinetics.
* An analogous Maple file `maple/functions.mpl`.
* A directory `results` that contains the following files:
    -  `investigated_models.csv` with all networks in [ODEbase](https://www.odebase.org/) (as of November 2, 2023) with at least one reaction.
    -  `nondegenerate_networks` with all networks from `investigated_models.csv` that admit a positive nondegenerate steady state.
    -  `degenerate_networks` with all networks from `investigated_models.csv` that have a positive steady states, but all of them are degenerate.
    -  `generic_local_acr` with all networks from `investigated_models.csv` that satisfy the following criteria:
       * admits nondegenerate positive steady states
       * is not of full rank (after removing nonparticipating species)
       * has generic local ACR in at least one speceis.

## Dependencies

The Julia portion of the code is based on Catalyst v14.4.1 and Oscar v1.1.1. For exact dependencies, see the file `Manifest.toml`.

The Maple portion of the code relies on Maple 2023.

## Julia example

We begin by loading the functions:

```julia
include("julia/functions.jl");
```

Consider the following isocitrate dehydrogenase that appears in Shinar–Feinberg's work on absolute concentration robustness, entered in catalyst format.

```julia
rn = @reaction_network begin 
    k1, X1 + X2 --> X3
    k2, X3 --> X1 + X2
    k3, X3 --> X1 + X4
    k4, X3 + X4 --> X5
    k5, X5 --> X3 + X4
    k6, X5 --> X2 + X3 
end;
```

The following command returns `true`, which means that the network admits positive steady states:

```julia
julia> is_consistent(rn)
true
```

The following command returns `true`, which means that there is a nondegenerate steady state with respect to its stoichiometric compatibility classes:

```julia
julia> has_nondegenerate_steady_state(rn, use_conservation_laws=true)
true
```

We check for generic local ACR with respect to the first and fourth species:

```julia
julia> generic_local_acr(rn, 1)
false

julia> generic_local_acr(rn, 4)
true
```

We could also do these checks on the level of the matrices that describe the associated augmented vertical system:

```julia
N = matrix(QQ, netstoichmat(rn))
B = matrix(ZZ, substoichmat(rn))
L = martix(QQ, conservationlaws(rn))

has_nondegenerate_zero(N, B, L)
generic_local_acr(N, B, 1)
generic_local_acr(N, B, 4)

```


## Maple example
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

The following command returns `false`; with the notation from equation (4.2) in the paper, this means that all zeros of $f_\kappa$ are degenerate for all $\kappa$:
```
HasNondegenerateZero(Gamma,B);
```
