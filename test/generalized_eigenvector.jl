# 1x1 matrix - trivial
alg_mult = [(1, 1)]
λ_basis, blocks = JordanForm.generalized_eigenvectors(O1 // 1, 1, 1)
expected_block_sizes = [(1, 1)]

@test length(blocks) == 1
@test expected_block_sizes == map(size, blocks)
λ_basis = map(v -> JordanForm.unwrap.(v), λ_basis)
@test λ_basis == [[1]] # FIXME: Test for subspace containment rather than exact eigenvectors

# 2x2 matrices
alg_mult = [(7 - JordanForm.symbolic_sqrt(41) * 1im, 1), (7 + JordanForm.symbolic_sqrt(41) * 1im, 1)]
λ_basis, blocks = JordanForm.generalized_eigenvectors(A1 // 1, alg_mult[1][1], alg_mult[1][2])
expected_block_sizes = [(1, 1)]

@test length(blocks) == 1
@test expected_block_sizes == map(size, blocks)
λ_basis = map(v -> JordanForm.unwrap.(v), λ_basis)
println(λ_basis[1])
@test λ_basis == [[1 // 6 + 41^(1 // 2) * 1im // 6, 1]] # FIXME: Test for subspace containment rather than exact eigenvectors

λ_basis, blocks = JordanForm.generalized_eigenvectors(A1 // 1, alg_mult[2][1], alg_mult[2][2])
expected_block_sizes = [(1, 1)]

@test length(blocks) == 1
@test expected_block_sizes == map(size, blocks)
λ_basis = map(v -> JordanForm.unwrap.(v), λ_basis)
@test λ_basis == [[1 // 6 - 41^(1 // 2) * 1im // 6, 1]] # FIXME: Test for subspace containment rather than exact eigenvectors
