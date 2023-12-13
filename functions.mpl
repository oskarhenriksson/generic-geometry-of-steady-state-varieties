# Returns a full-rank matrix W such that WN=0
ConservedQuantities := proc(N)

	local left_kernel, v:

	left_kernel := LinearAlgebra[NullSpace](LinearAlgebra[Transpose](N)):

	return Matrix([seq(convert(v,list),v in left_kernel)]):

end proc:

# Checks whether a matrix M has any nonzero p-by-p minor
HasNonzeroMinor := proc(M,p)
local chosen_rows, chosen_cols, minor;
	for chosen_rows in Iterator[Combination](LinearAlgebra[RowDimension](M),p) do
		for chosen_cols in Iterator[Combination](LinearAlgebra[ColumnDimension](M),p) do
			minor := expand(LinearAlgebra:-Determinant(M(chosen_rows+1,chosen_cols+1))):
			if minor <> 0 then
				return true:
			end if:
		end do:
	end do:
	return false:
end proc:

# Finds a matrix of full rank with the same row space as Gamma
RowReduce := proc(Gamma)
	local r:
	LinearAlgebra[RowSpace](Gamma):
	return <seq(<r>,r in %)>:
end proc:

# Checks whether Gamma has a positive vector in its kernel
# If Gamma is the stoichiometric matrix of a network, this corresponds to checking if the network is consistent
IsConsistent := proc(Gamma,max_time:=infinity)
	local existence_of_positive_flux, v, w, i, r;
	r := LinearAlgebra[ColumnDimension](Gamma):
	existence_of_positive_flux := {seq(v=0,v in convert(Gamma.Vector([seq(w[i],i=1..r)]),list)),seq(w[i]>0,i=1..r)};
	return SMTLIB[Satisfiable](existence_of_positive_flux ,logic="QF_NRA",timelimit=max_time);
end proc:

# Computes the matrix that appears in the rank condition (3.13) in the paper
# The rank of this matrix determines the existence of nondegenerate zeros of N.diag(k).x^B
DegeneracyMatrix := proc(N,B)
	local n,r,s,i,G;
	n := LinearAlgebra[RowDimension](B);
	r := LinearAlgebra[ColumnDimension](B);
	s := LinearAlgebra[Rank](N);
	LinearAlgebra[NullSpace](N):
	map(x->LinearAlgebra[Transpose](x),%):
	G := LinearAlgebra[Transpose](<op(%)>):
	return N.LinearAlgebra[DiagonalMatrix](G.Vector([seq(u[i],i=1..r-s)])).LinearAlgebra[Transpose](B):
end proc:

# Checks whether N.diag(k).x^B has a nondegenerate zero, where N has full rank and the same row space as Gamma
# If both this function and IsConsistent(Gamma) return true, then the network satisfies the equivalent conditions in Theorem 4.2
HasNondegenerateZero := proc(Gamma,B,number_of_attempts:=3,max_entry_size:=1000,max_time:=infinity)
	local n, r, s, G, degeneracy_matrix, N, i, j, ranks_at_random_steady_states:

	n := LinearAlgebra[RowDimension](B):
	r := LinearAlgebra[ColumnDimension](B):
	s := LinearAlgebra[Rank](Gamma):


	N := RowReduce(Gamma):

	degeneracy_matrix := DegeneracyMatrix(N,B):
	ranks_at_random_steady_states := [seq(LinearAlgebra[Rank](
					subs([seq(u[j]=rand(-max_entry_size..max_entry_size)(),j=1..r-s)],degeneracy_matrix)),
						i=1..number_of_attempts)]:
	
	if (s in ranks_at_random_steady_states) then
		return true:
	else
		if timelimit(max_time,HasNonzeroMinor(degeneracy_matrix,s)) then
			return true:
		else
			return false:
		end if:
	end if:
end proc:

# Checks whether a network with stoichiometric matrix Gamma and reactant matrix B has a nondegenerate steady state
# If both this function and IsConsistent(Gamma) return true, then the network satisfies the equivalent conditions in Theorem 4.5
HasNondegenerateSteadyState := proc(Gamma,B,number_of_attempts:=3,max_entry_size:=1000,max_time:=infinity)
	local n, r, s, N, G, W, degeneracy_matrix, augmented_degeneracy_matrix, i, j, ranks_at_random_steady_states:
	n := LinearAlgebra[RowDimension](B):
	r := LinearAlgebra[ColumnDimension](B):
	s := LinearAlgebra[Rank](Gamma):

	N := RowReduce(Gamma):

	degeneracy_matrix := DegeneracyMatrix(N,B):
	W := ConservedQuantities(Gamma):
	augmented_degeneracy_matrix := <degeneracy_matrix,W>:
	ranks_at_random_steady_states := [seq(LinearAlgebra[Rank](
		subs([seq(u[j]=rand(-max_entry_size..max_entry_size)(),j=1..r-s)],augmented_degeneracy_matrix)),
		i=1..number_of_attempts)]:
	if (n in ranks_at_random_steady_states) then
		return true:
	else
		if timelimit(max_time,HasNonzeroMinor(augmented_degeneracy_matrix,n)) then
			return true:
		else
			return false:
		end if:
	end if:
end proc:
