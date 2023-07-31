using JordanForm, SymbolicUtils
using Test

@testset "JordanForm.jl" begin
    include("matrix_data.jl")

    @testset include("charpoly.jl")
    @testset include("eigenvalue.jl")
    @testset include("generalized_eigenvector.jl")
    @testset include("jordan_form.jl")
end
