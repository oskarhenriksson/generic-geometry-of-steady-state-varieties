using Oscar
using Catalyst

# Basic functions

# Function to generate all p-by-p minors of matrix M
# (Will be included in future versions of Oscar!)
function minors_iterator(M::MatrixElem, k::Int)
    row_indices = AbstractAlgebra.combinations(1:nrows(M), k)
    col_indices = AbstractAlgebra.combinations(1:ncols(M), k)
    return (det(M[rows, cols]) for rows in row_indices for cols in col_indices)
end

function is_unit_vector(v::Vector)
    sum(v) == 1 && all(x -> x == 0 || x == 1, v)
end


# Vertically parametrized systems

# Checks the matrix of monomials for linearity
function vertically_parametrized_system(B::ZZMatrix, C::QQMatrix)
    Qk, k = rational_function_field(QQ, "k"=>1:size(B, 2))
    Qkx, x = polynomial_ring(Qk, "x"=>1:size(B, 1))
    return C*[k[i]*prod(x.^B[:,i]) for i in 1:ncols(B)]
end

function is_linear(B::ZZMatrix)
    for i in 1:ncols(B)
        col = B[:,i]
        if !is_unit_vector(col) && !is_zero(col)
            return false
        end
    end
    return true
end

# Checks whether N has a positive vector in its kernel
# If N is the stoichiometric matrix of a network, this corresponds to checking if the network is consistent
function nonempty_positive_kernel(N::Union{QQMatrix,ZZMatrix})
    inequalities = (-identity_matrix(QQ, ncols(N)), -ones(Int, ncols(N)))
    equalities = (N, zeros(Int, nrows(N)))
    P = polyhedron(inequalities, equalities)
    return is_feasible(P)
end


# Generic nondegeneracy checks

# Checks whether (C.diag(k).x^M, L*x-b) has a nondegenerate zero
function has_nondegenerate_zero(C::QQMatrix, M::ZZMatrix, L::QQMatrix=zero_matrix(QQ, 0, nrows(M));
    number_of_attempts::Int=3, max_entry_size::Int=1000, certify::Bool=true)
    C = (rref(C)[2])[1:rank(C), :]
    L = rref(L)[2]
    @req ncols(C) == ncols(M) "C and M need to have the same number of columns"
    @req ncols(L) == nrows(M) "L needs to have the same number of columns as M has rows"
    s = rank(C)
    G = kernel(C, side=:right)
    for _ in 1:number_of_attempts
        u = rand(-max_entry_size:max_entry_size, ncols(G))
        nrows(L) == 0 ? h = ones(Int, nrows(M)) : h = rand(-max_entry_size:max_entry_size, nrows(M))
        degeneracy_matrix = vcat(C * diagonal_matrix(G * u) * transpose(M) * diagonal_matrix(h), L)
        if rank(degeneracy_matrix) == s + nrows(L)
            return true
        end
    end
    if certify
        R, u, h = polynomial_ring(QQ, "u" => 1:ncols(G), "h" => 1:nrows(M))
        nrows(L) == 0 ? h = ones(Int, nrows(M)) : nothing
        symbolic_degeneracy_matrix = vcat(C * diagonal_matrix(G * u) * transpose(M) * diagonal_matrix(h), R.(L))
        return !all(is_zero, minors_iterator(symbolic_degeneracy_matrix, s + nrows(L)))
    else
        return false
    end
end

# ACR functions

function generic_local_acr(C::QQMatrix, M::ZZMatrix, i::Int; 
    number_of_attempts::Int=3, max_entry_size::Int=1000, certify::Bool=true)
    @req has_nondegenerate_zero(C, M) "The system needs to have a nondegenerate zero"
    return !has_nondegenerate_zero(C, M[setdiff(1:nrows(M), i), :], 
        number_of_attempts=number_of_attempts, max_entry_size=max_entry_size, certify=certify)
end


function local_acr_polynomial(C::QQMatrix, M::ZZMatrix, i::Int)
    steady_state_system = vertically_parametrized_system(M, C)
    R = parent(first(steady_state_system))
    x = gens(R)
    I = ideal(parent(first(steady_state_system)), steady_state_system)
    Isat = saturation(I, ideal(R, prod(x)))
    Ielim = eliminate(Isat, [x[j] for j in 1:nrows(M) if j != i])
    generators = gens(Ielim)
    @req length(generators) == 1 "The ideal should be principal"
    g = first(generators)
    return g
end


# CRNT functions

function is_consistent(rn::ReactionSystem)
    N = matrix(QQ, netstoichmat(rn))
    return nonempty_positive_kernel(N)
end

function has_nondegenerate_steady_state(rn::ReactionSystem; use_conservation_laws::Bool=false, 
    number_of_attempts::Int=3, max_entry_size::Int=1000, certify::Bool=true)
    N = matrix(QQ, netstoichmat(rn))
    B = matrix(ZZ, substoichmat(rn))
    if use_conservation_laws
        L = matrix(QQ, conservationlaws(rn))
    else
        L = zero_matrix(QQ, 0, nrows(B))
    end
    return has_nondegenerate_zero(N, B, L, number_of_attempts=number_of_attempts, max_entry_size=max_entry_size, certify=certify)
end

function generic_local_acr(rn::ReactionSystem, i::Int;
    number_of_attempts::Int=3, max_entry_size::Int=1000, certify::Bool=true)
    N = matrix(QQ, netstoichmat(rn))
    B = matrix(ZZ, substoichmat(rn))
    return generic_local_acr(N, B, i; number_of_attempts=number_of_attempts, max_entry_size=max_entry_size, certify=certify)
end

function local_acr_polynomial(rn::ReactionSystem, i::Int)
    N = matrix(QQ, netstoichmat(rn))
    B = matrix(ZZ, substoichmat(rn))
    return local_acr_polynomial(N, B, i)
end

has_linear_kinetics(rn::ReactionSystem) = is_linear(matrix(ZZ, substoichmat(rn)))
