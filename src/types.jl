export JordanBlock, JordanCanonicalForm, JordanFactorization
export eigenvalue, diaglength, blocks, nblocks

# Block
struct JordanBlock{T} <: AbstractMatrix{T}
    λ::T
    size::Integer
end
Base.eltype(::JordanBlock{T}) where {T} = T

eigenvalue(J::JordanBlock) = J.λ
diaglength(J::JordanBlock) = J.size

Base.size(J::JordanBlock) = (diaglength(J), diaglength(J))
Base.size(J::JordanBlock, d::Integer) = diaglength(J)

function Base.getindex(J::JordanBlock{T}, i::Integer, j::Integer) where {T}
    if i == j
        return J.λ
    elseif i == j - 1
        return one(T)
    else
        return zero(T)
    end
end

# Form
struct JordanCanonicalForm{T} <: AbstractMatrix{T}
    jordan_blocks::Vector{JordanBlock{T}}
end
JordanCanonicalForm(jordan_blocks::Vector{<:JordanBlock}) = JordanCanonicalForm{eltype(first(jordan_blocks))}(jordan_blocks)

Base.eltype(::JordanCanonicalForm{T}) where {T} = T
blocks(J::JordanCanonicalForm) = J.jordan_blocks
nblocks(J::JordanCanonicalForm) = length(blocks(J))

function Base.size(J::JordanCanonicalForm, d::Integer)
    n = sum(diaglength, blocks(J))
    return n
end
function Base.size(J::JordanCanonicalForm)
    n = size(J, 1)
    return (n, n)
end

function Base.getindex(J::JordanCanonicalForm{T}, i::Integer, j::Integer) where {T}
    # FIXME: Compute correct indexing.
    
    k = 1
    @inbounds for block in blocks(J)
        if k <= i < k + size(block, 1)
            if k <= j < k + size(block, 2)
                return block[i - k, j - k]
            else
                return zero(T)
            end
        else
            k += size(block, 1)
        end
    end

    throw(BoundsError(J, i))
end

# Factorization
struct JordanFactorization{S, T}
    S::AbstractMatrix{S}
    J::JordanCanonicalForm{T}
end

Base.iterate(F::JordanFactorization) = (F.S, Val(:form))
Base.iterate(F::JordanFactorization, ::Val{:form}) = (F.J, Val(:done))
Base.iterate(F::JordanFactorization, ::Val{:done}) = nothing
nblocks(F::JordanFactorization) = nblocks(F.J)