#  Generic consistency and nondegeneracy of vertically parametrized systems
This repository contains files for the manuscript "The generic geometry of steady state varieties"by Elisenda Feliu, Oskar Henriksson, and Beatriz Pascual-Escudero.

## File descriptions
The repository contains the following files:
* A Julia file `julia/functions.jl` that contains functions for testing whether a network admits positive nondegenerate steady states when modeled with (generalized) mass action kinetics.
* Two notebooks with examples:
  - `julia/IDH.ipynb` for the isocitrate dehydrogenase network in Example 4.1 of the paper.
  - `julia/167.ipynb` for the network `BIOMD0000000167` from BioModels/ODEbase discussed in Example 4.10 of the paper.
* A directory `results` that contains the following files:
    -  `investigated_models.csv` with all networks in [ODEbase](https://www.odebase.org/) (as of November 2, 2023) with at least one reaction.
    -  `nondegenerate_networks.csv` with all networks from `investigated_models.csv` that admit a positive nondegenerate steady state.
    -  `degenerate_networks.csv` with all networks from `investigated_models.csv` that have a positive steady states, but all of them are degenerate.
    -  `generic_local_acr.csv` with all networks from `investigated_models.csv` that satisfy the following criteria:
       * admits nondegenerate positive steady states
       * is not of full rank (after removing nonparticipating species)
       * has generic local ACR in at least one speceis.

## Dependencies

The Julia portion of the code is based on Catalyst v14.4.1 and Oscar v1.1.1. For exact dependencies, see the file `julia/Manifest.toml`.

The Maple portion of the code was written for Maple 2023.

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

```julia-repl
julia> is_consistent(rn)
true
```

The following command returns `true`, which means that there is a nondegenerate steady state with respect to its stoichiometric compatibility classes:

```julia-repl
julia> has_nondegenerate_steady_state(rn, use_conservation_laws=true)
true
```

We check for generic local ACR with respect to the first and fourth species:

```julia-repl
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
