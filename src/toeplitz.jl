struct NonSquareLowerTriangularToeplitz{T, V<:AbstractVector{T}} <: AbstractMatrix{T}
    v::V
    n_cols::Int64

    function NonSquareLowerTriangularToeplitz{T, V}(v::V, n_cols) where {T, V <: AbstractVector{T}}
        require_one_based_indexing(v)
        new{T, V}(v, n_cols)
    end
end

# TODO: Fix automatic n_cols
function NonSquareLowerTriangularToeplitz{T}(v::AbstractVector, n_cols) where T
    vT = convert(AbstractVector{T}, v)
    NonSquareLowerTriangularToeplitz{T, typeof(vT)}(vT, n_cols)
end
NonSquareLowerTriangularToeplitz(v::V, n_cols) where {T, V<:AbstractVector{T}} = NonSquareLowerTriangularToeplitz{T, V}(v, n_cols)

basetype(::NonSquareLowerTriangularToeplitz{T}) where {T <: NonSquareLowerTriangularToeplitz} = NonSquareLowerTriangularToeplitz

(==)(A::NonSquareLowerTriangularToeplitz, B::NonSquareLowerTriangularToeplitz) = A.n_cols == B.n_cols && A.v == B.v

function copyto!(A::NonSquareLowerTriangularToeplitz, B::NonSquareLowerTriangularToeplitz)
    copyto!(A.v, B.v)
    A.n_cols = B.n_cols
    A
end

parent(A::NonSquareLowerTriangularToeplitz) = A.v
basetype(x) = basetype(typeof(x))

function size(A::NonSquareLowerTriangularToeplitz)
    n = length(parent(A))
    (n, A.n_cols)
end

adjoint(A::NonSquareLowerTriangularToeplitz) = transpose(conj(A))
function zero!(A::NonSquareLowerTriangularToeplitz)
    fill!(parent(A), zero(eltype(A)))
    return A
end

function lmul!(x::Number, A::NonSquareLowerTriangularToeplitz)
    lmul!(x, parent(A))
    A
end
function rmul!(A::NonSquareLowerTriangularToeplitz, x::Number)
    rmul!(parent(A), x)
    A
end

iszero(A::NonSquareLowerTriangularToeplitz) = iszero(parent(A))

(*)(scalar::Number, C::NonSquareLowerTriangularToeplitz) = basetype(C)(scalar * parent(C))
(*)(C::NonSquareLowerTriangularToeplitz, scalar::Number) = basetype(C)(parent(C) * scalar)

AbstractMatrix{T}(A::NonSquareLowerTriangularToeplitz) where {T} = basetype(A){T}(AbstractVector{T}(A.v), A.n_cols)

for fun in (:zero, :conj, :copy, :-, :real, :imag)
    @eval $fun(A::NonSquareLowerTriangularToeplitz) = basetype(A)($fun(parent(A)), A.n_cols)
end

for op in (:+, :-)
    @eval function $op(A::NonSquareLowerTriangularToeplitz, B::NonSquareLowerTriangularToeplitz)
        n, m = length(A.v), length(B.v)
        if A.n_cols != B.n_cols || n != m
            throw("Incompatible matrix dimensions: ($n, $(A.n_cols)) and ($m, $(B.n_cols))")
        end

        return NonSquareLowerTriangularToeplitz($op(A.v, B.v), A.n_cols)
    end
end

# vc and vr
function getproperty(A::NonSquareLowerTriangularToeplitz, s::Symbol)
    if s == :vc
        getfield(A, :v)
    elseif s == :vr
        _firstnonzero(getfield(A, :v))
    elseif s == :uplo
        :L
    else
        getfield(A,s)
    end
end

function _firstnonzero(v::AbstractVector)
    w = zero(v)
    w[1] = v[1]
    w
end

# getindex
Base.@propagate_inbounds function getindex(A::NonSquareLowerTriangularToeplitz{T}, i::Integer, j::Integer) where T
    @boundscheck checkbounds(A, i, j)
    return i >= j ? A.v[i - j + 1] : zero(T)
end

function NonSquareLowerTriangularToeplitz{T}(A::AbstractMatrix) where T
    checksquare(A)
    return NonSquareLowerTriangularToeplitz{T}(_vc(A))
end
