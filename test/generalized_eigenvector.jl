# 1x1 matrix - trivial
# O1
alg_mult = [(1, 1)]
λ_basis, jordan_blocks = JordanForm.generalized_eigenvectors(O1, 1, 1)

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1.0]])

# 2x2 matrices
# A1
alg_mult = [
    (7 - JordanForm.symbolic_sqrt(41) * 1im, 1),
    (7 + JordanForm.symbolic_sqrt(41) * 1im, 1),
]

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(A1, alg_mult[1][1], alg_mult[1][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1 / 6 + 41^(1 / 2) * 1im / 6, 1]])

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(A1, alg_mult[2][1], alg_mult[2][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1 / 6 - 41^(1 / 2) * 1im / 6, 1]])

# A2
alg_mult = [(-5 - 6^(1 / 2), 1), (-5 + 6^(1 // 2), 1)]

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(A2, alg_mult[1][1], alg_mult[1][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[4 / 5 - 6^(1 / 2) / 5, 1]])

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(A2, alg_mult[2][1], alg_mult[2][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[4 / 5 + 6^(1 / 2) / 5, 1]])

# A3
alg_mult = [(1, 2)]

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(A3, alg_mult[1][1], alg_mult[1][2])

expected_block_sizes = [(1, 1), (1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1, 0], [0, 1]])
@test insubspace(λ_basis[2], [[1, 0], [0, 1]])
@test islinearlyindependent(λ_basis)

# A4
alg_mult = [(1, 2)]

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(A4, alg_mult[1][1], alg_mult[1][2])

expected_block_sizes = [(2, 2)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1, 0], [0, 1]])
@test insubspace(λ_basis[2], [[1, 0], [0, 1]])
@test islinearlyindependent(λ_basis)

# A5
alg_mult = [(-1, 1), (1, 1)]

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(A5, alg_mult[1][1], alg_mult[1][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[-1, 1]])

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(A5, alg_mult[2][1], alg_mult[2][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1, 1]])

# 3x3 matrices
# B1
alg_mult = [(1, 1), (2, 1), (3, 1)]

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(B1, alg_mult[1][1], alg_mult[1][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1, 0, 0]])

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(B1, alg_mult[2][1], alg_mult[2][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[0, 1, 0]])

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(B1, alg_mult[3][1], alg_mult[3][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[0, 0, 1]])

# B2
alg_mult = [(1, 1), (2, 1), (3, 1)]

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(B2, alg_mult[1][1], alg_mult[1][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1, 0, 0]])

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(B2, alg_mult[2][1], alg_mult[2][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1, 1, 0]])

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(B2, alg_mult[3][1], alg_mult[3][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1 / 2, 1, 1]])

# 4x4 matrices
# C1
alg_mult = [(1, 1), (2, 1), (3, 1), (4, 1)]

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(C1, alg_mult[1][1], alg_mult[1][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1, 0, 0, 0]])

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(C1, alg_mult[2][1], alg_mult[2][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[0, 1, 0, 0]])

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(C1, alg_mult[3][1], alg_mult[3][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[0, 0, 1, 0]])

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(C1, alg_mult[4][1], alg_mult[4][2])

expected_block_sizes = [(1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[0, 0, 0, 1]])

# C2
λ_basis, jordan_blocks = JordanForm.generalized_eigenvectors(C2, 2, 4)
expected_block_sizes = [(2, 2), (2, 2)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

# 5x5 matrices
# D1
alg_mult = [(-1, 2), (1, 3)]

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(D1, alg_mult[1][1], alg_mult[1][2])

expected_block_sizes = [(2, 2)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[-3, -2, 0, -1, 2], [1, 0, 0, 1, 0]])
@test insubspace(λ_basis[2], [[-3, -2, 0, -1, 2], [1, 0, 0, 1, 0]])
@test islinearlyindependent(λ_basis)

λ_basis, jordan_blocks =
    JordanForm.generalized_eigenvectors(D1, alg_mult[2][1], alg_mult[2][2])

expected_block_sizes = [(2, 2), (1, 1)]
@test length(jordan_blocks) == length(expected_block_sizes)
@test expected_block_sizes == map(size, jordan_blocks)

@test insubspace(λ_basis[1], [[1, 1, 1, 1, -1], [0, 1, 0, 1, 0]])
@test insubspace(λ_basis[2], [[1, 1, 1, 1, -1], [0, 1, 0, 1, 0]])
@test insubspace(λ_basis[3], [[-1, -1, 1, 0, 0]])
@test islinearlyindependent(λ_basis)
