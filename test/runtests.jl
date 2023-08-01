using JordanForm, SymbolicUtils
using Test

test_files = ["charpoly.jl", "eigenvalue.jl", "generalized_eigenvector.jl", "jordan_form.jl"]
include("matrix_data.jl")
    
for f = test_files
    @testset "$f" begin
        include(f)
    end
end
