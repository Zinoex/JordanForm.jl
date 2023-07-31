using JordanForm
using Test

test_files = ["eigenvalues_test.jl"]

for f = test_files
    @info "Test - $f"
    include(f)
end