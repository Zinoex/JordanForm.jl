p = JordanForm.charpoly(A1)
expected = [1, -14, 90]
@test iszero(expected - p)

p = JordanForm.charpoly(A2)
expected = [1, 10, 19]
@test iszero(expected - p)

p = JordanForm.charpoly(A3)
expected = [1, -2, 1]
@test iszero(expected - p)

p = JordanForm.charpoly(A4)
expected = [1, -2, 1]
@test iszero(expected - p)

p = JordanForm.charpoly(A5)
expected = [1, 0, -1]
@test iszero(expected - p)

p = JordanForm.charpoly(D1)
expected = [1, -1, -2, 2, 1, -1]
@test iszero(expected - p)