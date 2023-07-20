module JordanForm

using LinearAlgebra
using Graphs
using Symbolics

const IntOrRational = Union{Integer, Rational}

include("toeplitz.jl")
include("eigenvalues.jl")

export jordan_form

end
