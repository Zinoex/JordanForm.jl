# This is adapted from ToeplitzMatrices.jl with the modification that the Toeplitz matrix is not square.

import Base: size, getindex, parent

struct NonSquareLowerTriangularToeplitz{T, V<:AbstractVector{T}} <: AbstractMatrix{T}
    v::V
    n_cols::Int64

    function NonSquareLowerTriangularToeplitz{T, V}(v::V, n_cols) where {T, V <: AbstractVector{T}}
        Base.require_one_based_indexing(v)
        new{T, V}(v, n_cols)
    end
end

function NonSquareLowerTriangularToeplitz{T}(v::AbstractVector, n_cols) where T
    vT = convert(AbstractVector{T}, v)
    NonSquareLowerTriangularToeplitz{T, typeof(vT)}(vT, n_cols)
end
NonSquareLowerTriangularToeplitz(v::V, n_cols) where {T, V<:AbstractVector{T}} = NonSquareLowerTriangularToeplitz{T, V}(v, n_cols)

parent(A::NonSquareLowerTriangularToeplitz) = A.v
function size(A::NonSquareLowerTriangularToeplitz)
    n = length(parent(A))
    (n, A.n_cols)
end

# getindex
Base.@propagate_inbounds function getindex(A::NonSquareLowerTriangularToeplitz{T}, i::Integer, j::Integer) where T
    @boundscheck checkbounds(A, i, j)
    return i >= j ? A.v[i - j + 1] : zero(T)
end
