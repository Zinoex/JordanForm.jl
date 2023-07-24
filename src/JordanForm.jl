module JordanForm

using LinearAlgebra: checksquare, I, rank
using RowEchelon
using Graphs
using Symbolics

const IntOrRational = Union{Integer, Rational}

include("types.jl")
include("toeplitz.jl")
include("eigenvalues.jl")

export jordan_form

function jordan_form(M::AbstractMatrix{T}) where {T <: IntOrRational} 
    checksquare(M)
    
    eigs = radical_eigvals(M)
    eigs = algebraic_multiplicity(eigs)

    M = M // 1

    jordan_basis = []
    jordan_blocks = JordanBlock[]
    
    for (λ, alg_mul) in eigs
        E = M - λ * I

        # Exact nullspace computation
        λvecs = rref_nullspace(E)

        for v in eachcol(λvecs)
            push!(jordan_basis, v)
            s = 1

            for _ in 2:alg_mul
                v = rref_linsolve(E, v)[:, 1]

                push!(jordan_basis, v)
                s += 1

                if iszero(E * v)
                    break
                end
            end

            push!(jordan_blocks, JordanBlock(λ, s))
        end
    end
    
    basis_mat = reduce(hcat, jordan_basis)
    
    return JordanFactorization(basis_mat, JordanCanonicalForm(jordan_blocks))
end

rref_nullspace(A::AbstractMatrix{T}) where {T} = rref_linsolve(A, zeros(T, size(A, 1)))
function rref_linsolve(A::AbstractMatrix{T}, b::AbstractVector{T}) where {T}
    n, m = size(A)
    @assert n == m

    aug = hcat(A, b)
    aug, pivots = rref_with_pivots!(aug, 0)

    if last(pivots) == m + 1
        return nothing  # The system is inconsistent
    end

    frees = setdiff(collect(1:m), pivots)

    result = zeros(T, n, length(frees))
    for (k, j) in enumerate(frees)
        v = @view(result[:, k])
        v[j] = T(1)
        for (i, l) in enumerate(pivots)
            v[l] = aug[i, end]
            if l < j
                v[l] -= aug[i, j]
            end
        end
    end

    return result
end

end
