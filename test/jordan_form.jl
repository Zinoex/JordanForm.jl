# 1x1 matrix - trivial
# O1
F = jordan_form(O1)
S, J = F

@test size(S) == (1, 1)
@test size(J) == (1, 1)

expected_block_sizes = [(1, 1)]
@test nblocks(J) == length(expected_block_sizes)
@test nblocks(F) == length(expected_block_sizes)
@test expected_block_sizes == map(size, blocks(J))

expected_blocks = [(1, 1)]
@test expected_blocks == map(b -> (eigenvalue(b), (diaglength(b))), blocks(J))

@test islinearlyindependent(S)
@test O1 * S ≈ S * J

# 2x2 matrices
# A1
F = jordan_form(A1)
S, J = F

@test size(S) == (2, 2)
@test size(J) == (2, 2)

expected_block_sizes = [(1, 1), (1, 1)]
@test nblocks(J) == length(expected_block_sizes)
@test nblocks(F) == length(expected_block_sizes)
@test expected_block_sizes == map(size, blocks(J))

expected_blocks = [(7 - JordanForm.symbolic_sqrt(41) * 1im, 1), (7 + JordanForm.symbolic_sqrt(41) * 1im, 1)]
@test expected_blocks == map(b -> (eigenvalue(b), (diaglength(b))), blocks(J))

@test islinearlyindependent(S)
@test A1 * S ≈ S * J

# A2
F = jordan_form(A2)
S, J = F

@test size(S) == (2, 2)
@test size(J) == (2, 2)

expected_block_sizes = [(1, 1), (1, 1)]
@test nblocks(J) == length(expected_block_sizes)
@test nblocks(F) == length(expected_block_sizes)
@test expected_block_sizes == map(size, blocks(J))

expected_blocks = [(-5 - 6^(1/2), 1), (-5 + 6^(1//2), 1)]
@test expected_blocks == map(b -> (eigenvalue(b), (diaglength(b))), blocks(J))

@test islinearlyindependent(S)
@test A2 * S ≈ S * J

# A3
F = jordan_form(A3)
S, J = F

@test size(S) == (2, 2)
@test size(J) == (2, 2)

expected_block_sizes = [(1, 1), (1, 1)]
@test nblocks(J) == length(expected_block_sizes)
@test nblocks(F) == length(expected_block_sizes)
@test expected_block_sizes == map(size, blocks(J))

expected_blocks = [(1, 1), (1, 1)]
@test expected_blocks == map(b -> (eigenvalue(b), (diaglength(b))), blocks(J))

@test islinearlyindependent(S)
@test A3 * S ≈ S * J

# A4
F = jordan_form(A4)
S, J = F

@test size(S) == (2, 2)
@test size(J) == (2, 2)

expected_block_sizes = [(2, 2)]
@test nblocks(J) == length(expected_block_sizes)
@test nblocks(F) == length(expected_block_sizes)
@test expected_block_sizes == map(size, blocks(J))

expected_blocks = [(1, 2)]
@test expected_blocks == map(b -> (eigenvalue(b), (diaglength(b))), blocks(J))

@test islinearlyindependent(S)
@test A4 * S ≈ S * J

# A5
F = jordan_form(A5)
S, J = F

@test size(S) == (2, 2)
@test size(J) == (2, 2)

expected_block_sizes = [(1, 1), (1, 1)]
@test nblocks(J) == length(expected_block_sizes)
@test nblocks(F) == length(expected_block_sizes)
@test expected_block_sizes == map(size, blocks(J))

expected_blocks = [(-1, 1), (1, 1)]
@test expected_blocks == map(b -> (eigenvalue(b), (diaglength(b))), blocks(J))

@test islinearlyindependent(S)
@test A5 * S ≈ S * J

# 3x3 matrices
# B1
F = jordan_form(B1)
S, J = F

@test size(S) == (3, 3)
@test size(J) == (3, 3)

expected_block_sizes = [(1, 1), (1, 1), (1, 1)]
@test nblocks(J) == length(expected_block_sizes)
@test nblocks(F) == length(expected_block_sizes)
@test expected_block_sizes == map(size, blocks(J))

expected_blocks = [(1, 1), (2, 1), (3, 1)]
@test expected_blocks == map(b -> (eigenvalue(b), (diaglength(b))), blocks(J))

@test islinearlyindependent(S)
@test B1 * S ≈ S * J

# B2
F = jordan_form(B2)
S, J = F

@test size(S) == (3, 3)
@test size(J) == (3, 3)

expected_block_sizes = [(1, 1), (1, 1), (1, 1)]
@test nblocks(J) == length(expected_block_sizes)
@test nblocks(F) == length(expected_block_sizes)
@test expected_block_sizes == map(size, blocks(J))

expected_blocks = [(1, 1), (2, 1), (3, 1)]
@test expected_blocks == map(b -> (eigenvalue(b), (diaglength(b))), blocks(J))

@test islinearlyindependent(S)
@test B2 * S ≈ S * J

# 4x4 matrices
# C1
F = jordan_form(C1)
S, J = F

@test size(S) == (4, 4)
@test size(J) == (4, 4)

expected_block_sizes = [(1, 1), (1, 1), (1, 1), (1, 1)]
@test nblocks(J) == length(expected_block_sizes)
@test nblocks(F) == length(expected_block_sizes)
@test expected_block_sizes == map(size, blocks(J))

expected_blocks = [(1, 1), (2, 1), (3, 1), (4, 1)]
@test expected_blocks == map(b -> (eigenvalue(b), (diaglength(b))), blocks(J))

@test islinearlyindependent(S)
@test C1 * S ≈ S * J

# 5x5 matrices
# D1
F = jordan_form(D1)
S, J = F

@test size(S) == (5, 5)
@test size(J) == (5, 5)

expected_block_sizes = [(2, 2), (2, 2), (1, 1)]
@test nblocks(J) == length(expected_block_sizes)
@test nblocks(F) == length(expected_block_sizes)
@test expected_block_sizes == map(size, blocks(J))

expected_blocks = [(-1, 2), (1, 2), (1, 1)]
@test expected_blocks == map(b -> (eigenvalue(b), (diaglength(b))), blocks(J))

@test islinearlyindependent(S)
@test D1 * S ≈ S * J