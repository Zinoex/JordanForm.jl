# JordanForm

[![Build Status](https://github.com/Zinoex/JordanForm.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Zinoex/JordanForm.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/Zinoex/JordanForm.jl/graph/badge.svg?token=ENG7LBLR1J)](https://codecov.io/gh/Zinoex/JordanForm.jl)

[JordanForm.jl](https://github.com/Zinoex/JordanForm.jl) is an implementation for computing the Jordan form, also known as Jordan normal form and Jordan canonical form, __intended for educational use__. The Jordan form $J$ is a useful almost diagonal structure that is always computable, even if a matrix $A$ is not diagonaliziable. It is commonly used in the study of linear algebra and differential equations. In particular, if we consider a linear system of differential equations $\dot{x} = Ax$, then the Jordan form of $A$ can be used to find the general solution $x(t) = e^{At}x_0 = Se^{Jt}S^{-1}x_0$ of the system where $S$ is the transformation matrix between $A$ and $J$. That is, $AS = SJ$. Furthermore, there exists an analytic formula for $e^{Jt}$. The transformation matrix $S$ is computed as the concatenation of the generalized eigenvectors of $A$.

However, the structure of $J$ is very unstable under perturbation, as it depends on repeated eigenvalues of $A$. Therefore, it is not recommended to use the Jordan form for numerical computations. So why did we write a library to compute the Jordan form? The Jordan form is a useful tool for studying and understanding matrices and linear systems - and it is fun to write code! 

You can read more about the Jordan form on [Wikipedia](https://en.wikipedia.org/wiki/Jordan_normal_form).

> [!WARNING]
> The Jordan form is not numerically stable. Attempting to use this with floating point numbers may yield unexpected results. This package is intended for educational use only.

## Usage
As previously noted, the Jordan form is not numerically stable. Therefore, we restrict the input to integer and rational matrices. The Jordan form of a matrix $A$ is computed using the function `jordan_form(A)`. The function returns a factorization struct `F = (S, J)` where `J` is the Jordan form of `A` and `S` is the transformation matrix such that `A * S = S * J`.

```julia
using JordanForm 

A = [
    1 1 0;
    0 2 1;
    0 0 3
]

F = jordan_form(A)  # or
S, J = jordan_form(A)

@assert A * S ≈ S * J
```