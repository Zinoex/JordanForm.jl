function radical_eigenvalues(A::AbstractMatrix{T}) where {T <: IntOrRational}
    n = checksquare(A)
    A = convert.(Rational{eltype(A)}, A)

    A, λs = cofactoreigs(A)
    
    p = charpoly(A)
end

function charpoly(A::AbstractMatrix{<:Rational})
    n = size(A, 1)

    @syms λ::Real   # Assume real eigenvalues

    return det(A - Diagonal(fill(λ, n)))
end

function reduce_lower_order_terms()

end

function cofactoreigs(A::AbstractMatrix{<:Rational})
    d = axes(A, 1)

    λ_ids = filter(i -> iscofactoreig(A, i), d)
    λs = A[λ_ids, λ_ids]
    A = A[Not(λ_ids), Not(λ_ids)]

    return A, λs
end

iscofactoreig(A, i) = almost_zero_row(A, i) || almost_zero_col(A, i)
almost_zero_row(A, i) = all(iszero, A[i, Not(i)])
almost_zero_col(A, i) = all(iszero, A[Not(i), i])