
λ = JordanForm.radical_eigvals(A1)
expected = [7 + 41 // 2 * 1im, 7 - 41 // 2 * 1im]
println(λ)
@test iszero(simplify(expected - λ))

# λ = JordanForm.radical_eigvals(A2)
# expected = [1, 10, 19]
# @test iszero(expected - p)

# λ = JordanForm.radical_eigvals(A3)
# expected = [1, -2, 1]
# @test iszero(expected - p)

# λ = JordanForm.radical_eigvals(A4)
# expected = [1, -2, 1]
# @test iszero(expected - p)

# λ = JordanForm.radical_eigvals(A5)
# expected = [1, 0, -1]
# @test iszero(expected - p)

# λ = JordanForm.radical_eigvals(D1)
# expected = [1, -1, -2, 2, 1, -1]
# @test iszero(expected - p)