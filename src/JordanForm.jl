module JordanForm

using LinearAlgebra, BlockDiagonals, RowEchelon, LinearAlgebraX
using Graphs
using Symbolics

const IntOrRational = Union{Integer, Rational}

include("toeplitz.jl")
include("eigenvalues.jl")

export JordanBlock, JordanCanonicalForm, JordanFactorization
export jordan_form

struct JordanBlock{T} <: AbstractMatrix{T}
    λ::T
    size::Integer
end
Base.size(J::JordanBlock) = (J.size, J.size)
Base.size(J::JordanBlock, d::Integer) = (J.size, J.size)
function Base.getindex(J::JordanBlock, i::Integer, j::Integer)
    if i == j
        return J.λ
    elseif i == j + 1
        return one(T)
    else
        return zero(T)
    end
end

struct JordanCanonicalForm{T} <: AbstractMatrix{T}
    jordan_blocks::Vector{JordanBlock{T}}
end
function Base.size(J::JordanCanonicalForm, d::Integer)
    n = sum(b -> size(b, 1), J.jordan_blocks)
    return n
end
function Base.size(J::JordanCanonicalForm)
    n = size(J, 1)
    return (n, n)
end
function Base.getindex(J::JordanCanonicalForm, i::Integer, j::Integer)
    # FIXME: Compute correct indexing.
    block = J.jordan_blocks[findfirst(b -> i <= size(b, 1), J.jordan_blocks)]
    return block[i, j]
end

mutable struct JordanFactorization{T}
    basis::AbstractMatrix{T}
    form::JordanCanonicalForm{T}
end

Base.iterate(J::JordanFactorization) = (J.basis, Val(:form))
Base.iterate(J::JordanFactorization, ::Val{:form}) = (J.form, Val(:done))
Base.iterate(J::JordanFactorization, ::Val{:done}) = nothing

jordan_form(M; kwargs...) = trytoint.(_jordan_form(M; kwargs...))
function _jordan_form(M::AbstractMatrix{T}) where {T <: IntOrRational} 
    LinearAlgebra.checksquare(M)
    
    eigs = radical_eigvals(M)
    eigs = algebraic_multiplicity(eigs)
            
    # f((eig, size)) = jordan_block(eig, size)
    # blocks = map(f, block_structure)
    # jordan_mat = BlockDiagonal(blocks)

    
    jordan_basis = []
    jordan_blocks = []
    
    for (λ, alg_mul) in eigs
        λ_basis = []

        println(λ, " ", alg_mul)
        display(nullspacex(M - λ * I))

        eigvecs = nullspacex(M - λ * I)
        new_vecs = []

        for v in eachcol(eigvecs)

            vec = eig_mat(λ, i) * [:, 1]
            push!(λ_basis, vec)
            push!(jordan_basis, vec)
        end
    end
    
    basis_mat = hcat(jordan_basis...)
    
    return JordanFactorization(basis_mat, JordanCanonicalForm(jordan_blocks))
end

end
