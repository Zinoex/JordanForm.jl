using JordanForm, SymbolicUtils, LinearAlgebra
using Test

function insubspace(v, subspace)
    subspace = reduce(hcat, subspace)
    return rank(hcat(subspace, v)) == rank(subspace)
end

function islinearlyindependent(vectors::AbstractVector{<:AbstractVector})
    return rank(reduce(hcat, vectors)) == length(vectors)
end

function islinearlyindependent(vectors::AbstractMatrix)
    return rank(vectors) == size(vectors, 2)
end

include("matrix_data.jl")
    
test_files = ["charpoly.jl", "eigenvalue.jl", "generalized_eigenvector.jl", "jordan_form.jl"]
for f = test_files
    @testset "$f" begin
        include(f)
    end
end
