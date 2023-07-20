module JordanForm

using LinearAlgebra, InvertedIndices
using Symbolics

const IntOrRational = Union{Integer, Rational}

include("eigenvalues.jl")

export jordan_form

end
