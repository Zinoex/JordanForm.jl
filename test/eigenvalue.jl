# 2x2 matrices
λ = JordanForm.radical_eigvals(A1) 
λ = simplify.(λ)

expected = [7 - 41^(1 // 2) * 1im, 7 + 41^(1 // 2) * 1im]
expected_mult = [(e, 1) for e in expected]
@test isequal(expected, λ)
@test isequal(JordanForm.algebraic_multiplicity(λ), expected_mult)

λ = JordanForm.radical_eigvals(A2)
λ = simplify.(λ)
expected = [-5 - 6^(1//2), -5 + 6^(1//2)]
expected_mult = [(e, 1) for e in expected]
@test isequal(expected, λ)
@test isequal(JordanForm.algebraic_multiplicity(λ), expected_mult)

λ = JordanForm.radical_eigvals(A3)
expected = [1, 1]
expected_mult = [(1, 2)]
@test isequal(expected, λ)
@test isequal(JordanForm.algebraic_multiplicity(λ), expected_mult)

λ = JordanForm.radical_eigvals(A4)
expected = [1, 1]
expected_mult = [(1, 2)]
@test isequal(expected, λ)
@test isequal(JordanForm.algebraic_multiplicity(λ), expected_mult)

λ = JordanForm.radical_eigvals(A5)
λ = simplify.(λ)
expected = [-1, 1]
expected_mult = [(-1, 1), (1, 1)] 
@test isequal(expected, λ)
@test isequal(JordanForm.algebraic_multiplicity(λ), expected_mult)

# # 3x3 matrices
λ = JordanForm.radical_eigvals(B1)
expected = [1, 2, 3]
expected_mult = [(e, 1) for e in expected]
@test isequal(expected, λ)
@test isequal(JordanForm.algebraic_multiplicity(λ), expected_mult)

λ = JordanForm.radical_eigvals(B2)
expected = [1, 2, 3]
expected_mult = [(e, 1) for e in expected]
@test isequal(expected, λ)
@test isequal(JordanForm.algebraic_multiplicity(λ), expected_mult)

# # 4x4 matrices
λ = JordanForm.radical_eigvals(C1)
expected = [1, 2, 3, 4]
expected_mult = [(e, 1) for e in expected]
@test isequal(expected, λ)
@test isequal(JordanForm.algebraic_multiplicity(λ), expected_mult)

# # 5x5 matrices
λ = JordanForm.radical_eigvals(D1)
λ = simplify.(λ)
expected = [-1, -1, 1, 1, 1]
expected_mult = [(-1, 2), (1, 3)]
@test isequal(expected, λ)
@test isequal(JordanForm.algebraic_multiplicity(λ), expected_mult)
