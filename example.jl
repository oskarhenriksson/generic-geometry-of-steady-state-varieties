using Oscar
using Catalyst

include("functions.jl")

rn = @reaction_network begin 
    k1, X1 + X2 --> X3
    k2, X3 --> X1 + X2
    k3, X3 --> X1 + X4
    k4, X3 + X4 --> X5
    k5, X5 --> X3 + X4
    k6, X5 --> X2 + X3 
end;

N = matrix(QQ,netstoichmat(rn))
B = matrix(ZZ,substoichmat(rn))

has_nondegenerate_zero(N,B)
generic_local_acr(N,B,1)
generic_local_acr(N,B,2)
generic_local_acr(N,B,3)
generic_local_acr(N,B,4)
generic_local_acr(N,B,5)

is_consistent(rn)
has_nondegenerate_steady_state(rn, use_conservation_laws=true)
generic_local_acr(rn, 1)
generic_local_acr(rn, 2)
generic_local_acr(rn, 3)
generic_local_acr(rn, 4)
generic_local_acr(rn, 5)