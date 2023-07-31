module JordanForm

using LinearAlgebra: checksquare, I, rank
# using RowEchelon
using Graphs
using Symbolics

const IntOrRational = Union{Integer, Rational}

export radical_eigvals, algebraic_multiplicity

include("types.jl")
include("toeplitz.jl")
include("eigenvalues.jl")

export jordan_form

function jordan_form(M::AbstractMatrix{T}) where {T <: IntOrRational} 
    F = _jordan_form(M)

    # TODO: Convert back into most specific format
    basis = trytoint.(F.basis)
    F = JordanFactorization(basis, F.form)

    return F
end
trytoint(x::Rational) = isinteger(x) ? Int64(x) : x
trytoint(x) = x

function _jordan_form(M::AbstractMatrix{T}) where {T <: IntOrRational}     
    eigs = radical_eigvals(M)
    eigs = algebraic_multiplicity(eigs)

    M = M // 1

    jordan_basis = []
    jordan_blocks = JordanBlock[]
    
    for (λ, alg_mul) in eigs
        basis, blocks = generalized_eigenvectors(M, λ, alg_mul)
        append!(jordan_basis, basis)
        append!(jordan_blocks, blocks)
    end

    # return nothing
    
    jordan_basis = reduce(hcat, jordan_basis)
    
    return JordanFactorization(jordan_basis, JordanCanonicalForm(jordan_blocks))
end

function generalized_eigenvectors(M, λ, alg_mul)
    n = checksquare(M)

    E = M - λ * I

    # Compute block structure
    powers = Vector{typeof(E)}(undef, alg_mul)
    powers[1] = E
    nilpotent_power = 1

    for i in 2:alg_mul
        powers[i] = E * powers[i - 1]
        nilpotent_power = i

        if rank(powers[i]) == 0
            break
        end
    end
    resize!(powers, nilpotent_power)

    nullity = map(p -> n - rank(p), powers)
    blocks = blocks_from_nullity(nullity)

    # Compute λ basis
    λ_basis = []
    for size in blocks
        if size == 1
            v = rref_nullspace(powers[size])[:, 1]
            push!(λ_basis, v)
        else 
            null_big = rref_nullspace(powers[size])
            null_small = rref_nullspace(powers[size - 1])
            if !isempty(λ_basis)
                null_small = reduce(hcat, [null_small; λ_basis])
            end
    
            v = pick_vec(null_big, null_small)
            vs = Vector{AbstractVector{eltype(v)}}(undef, size)
            vs[1] = v

            for i in 2:size
                v = E * v
                vs[i] = v
            end

            reverse!(vs)
            append!(λ_basis, vs)
        end
    end

    blocks = JordanBlock.(λ, blocks)

    λ_basis, blocks
end

function blocks_from_nullity(nullity)
    nullity = [0; nullity; nullity[end]]

    prev, cur, next = nullity[1:end - 2], nullity[2:end - 1], nullity[3:end]
    block_sizes = @. 2cur - prev - next

    blocks = [i for (i, size_num) in enumerate(block_sizes) for _ in 1:size_num]
    reverse!(blocks)

    return blocks
end

rref_nullspace(A::AbstractMatrix{T}) where {T} = rref_linsolve(A, zeros(T, size(A, 1)))
function rref_linsolve(A::AbstractMatrix{T}, b::AbstractVector{T}) where {T}
    n, m = size(A)
    @assert n == m

    aug = hcat(A, b)
    aug, pivots = rref_with_pivots!(aug)

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

function pick_vec(null_big, null_small)
    @assert size(null_small, 2) >= 1
    n = size(null_small, 2) + 1

    for v in eachcol(null_big)
        if rank(hcat(null_small, v)) == n
            return v
        end
    end
end


# function scale_eigenvector(v::AbstractVector{<:Rational})
#     factor = lcm(denominator.(v))
#     v = v * factor

#     return v
# end
# scale_eigenvector(v::AbstractVector) = v

### Attributed to RowEchelon.jl
### Modified to use rational division.
function rref_with_pivots!(A::Matrix{T}) where T
    nr, nc = size(A)
    pivots = Vector{Int64}()
    i = j = 1
    while i <= nr && j <= nc
        (m, mi) = findmax(abs.(A[i:nr,j]))
        mi = mi+i - 1
        if m <= 0
            if 0 > 0
                A[i:nr,j] .= zero(T)
            end
            j += 1
        else
            for k=j:nc
                A[i, k], A[mi, k] = A[mi, k], A[i, k]
            end
            d = A[i,j]
            for k = j:nc
                A[i,k] //= d
            end
            for k = 1:nr
                if k != i
                    d = A[k,j]
                    for l = j:nc
                        A[k,l] -= d*A[i,l]
                    end
                end
            end
            append!(pivots,j)
            i += 1
            j += 1
        end
    end
    return A, pivots
end
rref_with_pivots(A::Matrix{T}) where {T} = rref_with_pivots!(copy(A))

end
