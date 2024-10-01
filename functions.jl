using Oscar
using Catalyst

# Function to generate all p-by-p minors of matrix M
# (Will probably be included in latest version of Oscar!)
function minors_iterator(M::MatrixElem, k::Int)
    row_indices = AbstractAlgebra.combinations(1:nrows(M), k)
    col_indices = AbstractAlgebra.combinations(1:ncols(M), k)
    return ( det(M[rows,cols]) for rows in row_indices for cols in col_indices )
end

# Checks whether N has a positive vector in its kernel
# If N is the stoichiometric matrix of a network, this corresponds to checking if the network is consistent
function nonempty_positive_kernel(N::Union{QQMatrix,ZZMatrix})
    inequalities = (-identity_matrix(QQ, ncols(N)), -ones(Int,ncols(N)))
    equalities = (N, zeros(Int, nrows(N)))
    P = polyhedron(inequalities, equalities)
    return is_feasible(P)
end

# Checks whether (C.diag(k).x^M, L*x-b) has a nondegenerate zero
function has_nondegenerate_zero(C::QQMatrix, M::ZZMatrix, L::QQMatrix=zero_matrix(QQ,0,nrows(M)); 
        number_of_attempts::Int=3, max_entry_size::Int=1000)
    C = (rref(C)[2])[1:rank(C),:]
    L = rref(L)[2]
    @req ncols(C) == ncols(M) "C and M need to have the same number of columns"
    @req ncols(L) == nrows(M) "L needs to have the same number of columns as M has rows"
    s = rank(C)
    G = kernel(C, side=:right)
    for _ in 1:number_of_attempts
        u = rand(-max_entry_size:max_entry_size, ncols(G))
        nrows(L) == 0 ? h = ones(Int, nrows(M)) : h = rand(-max_entry_size:max_entry_size, nrows(M))
        degeneracy_matrix = vcat(C*diagonal_matrix(G*u)*transpose(M)*diagonal_matrix(h), L)
        if rank(degeneracy_matrix) == s + nrows(L)
            return true
        end
    end
    R, u, h = polynomial_ring(QQ, "u"=>1:ncols(G), "h"=>1:nrows(M))
    nrows(L) == 0 ? h = ones(Int, nrows(M)) : nothing
    symbolic_degeneracy_matrix = vcat(C*diagonal_matrix(G*u)*transpose(M)*diagonal_matrix(h), R.(L))
    return !all(is_zero, minors_iterator(symbolic_degeneracy_matrix, s+nrows(L)))
end

function has_nondegenerate_zero(rn::ReactionSystem; use_conservation_laws::Bool = false)
    N = matrix(QQ,netstoichmat(rn))
    B = matrix(ZZ,substoichmat(rn))
    if use_conservation_laws
        L = matrix(QQ,conservationlaws(rn))
    else
        L = zero_matrix(QQ,0,nrows(B))
    end
    return has_nondegenerate_zero(N,B,L)
end


function generic_local_acr(C::QQMatrix, M::ZZMatrix, i::Int)
    @req has_nondegenerate_zero(C,M) "The system needs to have a nondegenerate zero"
    return !has_nondegenerate_zero(C, M[setdiff(1:nrows(M),i),:])
end

function generic_local_acr(rn::ReactionSystem, i::Int)
    N = matrix(QQ,netstoichmat(rn))
    B = matrix(ZZ,substoichmat(rn))
    return generic_local_acr(N,B,i)
end
