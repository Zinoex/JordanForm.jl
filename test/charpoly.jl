# 1x1 matrix - trivial
p = JordanForm.charpoly(O1)
expected = [1, -1]
@test expected == p

# 2x2 matrices
p = JordanForm.charpoly(A1)
expected = [1, -14, 90]
@test expected == p

p = JordanForm.charpoly(A2)
expected = [1, 10, 19]
@test expected == p

p = JordanForm.charpoly(A3)
expected = [1, -2, 1]
@test expected == p

p = JordanForm.charpoly(A4)
expected = [1, -2, 1]
@test expected == p

p = JordanForm.charpoly(A5)
expected = [1, 0, -1]
@test expected == p

# 3x3 matrices
p = JordanForm.charpoly(B1)
expected = [1, -6, 11, -6]
@test expected == p

p = JordanForm.charpoly(B2)
expected = [1, -6, 11, -6]
@test expected == p

# 4x4 matrices
p = JordanForm.charpoly(C1)
expected = [1, -10, 35, -50, 24]
@test expected == p

# 5x5 matrices
p = JordanForm.charpoly(D1)
expected = [1, -1, -2, 2, 1, -1]
@test expected == p

# Particular polynomials
p = [1, -3, 3, -1]
@test JordanForm.iscubic(p)
@test isequal(JordanForm.simplify.(JordanForm.cubicroots(p)), [1 // 1, 1 // 1, 1 // 1])

p = [1, -4, 6, -4, 1]
@test JordanForm.isquartic(p)
@test isequal(
    JordanForm.simplify.(JordanForm.quarticroots(p)),
    [1 // 1, 1 // 1, 1 // 1, 1 // 1],
)

p = [1, 2, 3, 4, 5, 6, 7]
try
    JordanForm.charpoly_roots(p)
catch e
    @test e == "Unable to find symbolic roots of polynomial"
end
