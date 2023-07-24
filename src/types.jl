export JordanBlock, JordanCanonicalForm, JordanFactorization

# Block
struct JordanBlock{T} <: AbstractMatrix{T}
    λ::T
    size::Integer
end
Base.eltype(::JordanBlock{T}) where {T} = T

Base.size(J::JordanBlock) = (J.size, J.size)
Base.size(J::JordanBlock, d::Integer) = J.size

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

function Base.size(J::JordanCanonicalForm, d::Integer)
    n = sum(b -> size(b, 1), J.jordan_blocks)
    return n
end
function Base.size(J::JordanCanonicalForm)
    n = size(J, 1)
    return (n, n)
end

function Base.getindex(J::JordanCanonicalForm{T}, i::Integer, j::Integer) where {T}
    # FIXME: Compute correct indexing.
    
    k = 1
    @inbounds for block in J.jordan_blocks
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
    basis::AbstractMatrix{S}
    form::JordanCanonicalForm{T}
end

Base.iterate(J::JordanFactorization) = (J.basis, Val(:form))
Base.iterate(J::JordanFactorization, ::Val{:form}) = (J.form, Val(:done))
Base.iterate(J::JordanFactorization, ::Val{:done}) = nothing
