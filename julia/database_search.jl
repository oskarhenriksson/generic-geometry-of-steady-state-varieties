include("functions.jl")

# Set parameters
max_species = Inf

# Define the directory containing the models
#cd(@__DIR__)
const model_directory = "../odebase_2023"

# Retrieve the list of models in the directory
list_of_models = readdir(model_directory)

number_of_models = length(list_of_models)

# Initialize lists to categorize the models
degenerate = String[]
nondegenerate_wrt_S = String[]
nondegenerate_not_wrt_S = String[]
inconsistent = String[]
missing_files = String[]
error = String[]
gen_local_acr = String[]
has_irrelevant_species = String[]
full_rank = String[]
linear = String[]

function write_both(file, msg)
    println(msg)         # Write to terminal
    println(file, msg)   # Write to file
end

# Function for reading files of the ODE base format
read_matrix = function (path)
    file_content = read(path, String)
    cleaned_content = replace(file_content, r"[<>\;]" => "")
    lines = split(cleaned_content, "\n")
    return [parse.(Int, split(line, ",")) for line in lines if !isempty(line)]
end

# Function for collecting the kinetic and stoichiometric matrices given a model id
function collect_matrices(model_id)
    B_path = joinpath(model_directory, model_id, "kinetic_matrix.txt")
    N_path = joinpath(model_directory, model_id, "reconfigured_stoichiometric_matrix.txt")
    B = matrix(ZZ, read_matrix(B_path))
    N = matrix(QQ, read_matrix(N_path))
    return B, N
end

# Iterate over each model
open("output.txt", "w") do file

    for model_id in list_of_models
        write_both(file,"")
        write_both(file, model_id)

        # Read and parse the kinetic matrix (B) and stoichiometric matrix (N)
        B = nothing
        N = nothing
        try
            B, N = collect_matrices(model_id)
        catch e
            write_both(file, "Error reading matrices for model: $model_id")
            println(e)
            push!(error, model_id)
            continue  # Skip to the next model
        end

        # Check for irrelevant species
        irrelevant_species = [i for i=1:nrows(B) if (all(iszero, B[i,:]) && all(iszero, N[i,:]))]
        if !isempty(irrelevant_species)
            write_both(file, "Irrelevant species (ignored in subsequent analysis): $irrelevant_species")
            push!(has_irrelevant_species, model_id)
            continue
        end

        # Ignore irrelevant species
        N = N[setdiff(1:nrows(N), irrelevant_species), :]
        B = B[setdiff(1:nrows(B), irrelevant_species), :]


        # Check for linear kinetics
        if has_linear_kinetics(B)
            write_both(file, "Linear kinetics")
            push!(linear, model_id)
        end

        # Determine dimensions and ranks
        n = size(B, 1)
        r = size(B, 2)

        s = rank(N)
        d = n - s

        if d == 0
            push!(full_rank, model_id)
            write_both(file, "(Effectively) full rank")
            continue
        end

        # Skip models exceeding the maximum number of species
        if n > max_species
            continue
        end

        if !nonempty_positive_kernel(N)
            write_both(file, "Inconsistent")
            push!(inconsistent, model_id)
            continue  
        else
            write_both(file, "Consistent")
        end

        L = kernel(N, side=:left)

        nondeg_flag = false

        if has_nondegenerate_zero(N, B, L)
            push!(nondegenerate_wrt_S, model_id)
            write_both(file, "Nondegenerate wrt stoichiometric subspace")
            nondeg_flag = true
        elseif has_nondegenerate_zero(N, B)
            push!(nondegenerate_not_wrt_S, model_id)
            write_both(file, "Nondegenerate (but not wrt stoichiometric subspace)")
            nondeg_flag = true
        else
            push!(degenerate, model_id)
            write_both(file, "Degenerate")
        end

        if nondeg_flag && rank(N) < n
            if n < 2
                continue
            end
            gen_local_acr_counter = 0
            for i = 1:n
                if generic_local_acr(N, B, i, certify=false, number_of_attempts=5)
                    write_both(file, "Generic local ACR wrt X$i")
                    gen_local_acr_counter = gen_local_acr_counter + 1
                end
            end
            if gen_local_acr_counter > 0
                push!(gen_local_acr, model_id)
            end
        end
    end
end

# Summerize the results
println("\nDegenerate networks:\n", degenerate)
println(length(degenerate))

println("\nNondegenerate networks wrt S:\n", nondegenerate_wrt_S)
println(length(nondegenerate_wrt_S))

println("\nNondegenerate networks, but not wrt S:\n", nondegenerate_not_wrt_S)
println(length(nondegenerate_not_wrt_S))

println("\nInconsistent networks:\n", inconsistent)
println(length(inconsistent))

println("\nNetworks with missing files:\n", missing_files)
println(length(missing_files))

println("\nNetworks with errors:\n", error)
println(length(error))

println("\nNetworks with generic local ACR:\n", gen_local_acr)
println(length(gen_local_acr))

println("\nNetworks with irrelevant species:\n", has_irrelevant_species)
println(length(has_irrelevant_species))

println("\nNetworks with full rank:\n", full_rank)
println(length(full_rank))

println("\nNetworks with linear kinetics:\n", linear)
println(length(linear))