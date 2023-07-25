function radical_eigvals(A::AbstractMatrix{T}) where {T <: IntOrRational}
    checksquare(A)
    
    g = DiGraph{UInt16}(A)
    blocks = strongly_connected_components(g)

    eigs = mapreduce(vcat, blocks) do block
        # Copy the block into a new matrix for efficiency for computing the characteristic polynomial
        block_matrix = A[block, block]

        # Find the eigenvalues of the block
        p = charpoly(block_matrix)
        eigs = charpoly_roots(p)

        return eigs
    end

    sort!(eigs, by=λ -> (real(λ), imag(λ)))

    return eigs
end

function algebraic_multiplicity(eigs::AbstractVector{T}) where T
    mult_eigs = Tuple{T, Int64}[]

    @inbounds for λ in eigs
        if isempty(mult_eigs) || mult_eigs[end][1] != λ
            push!(mult_eigs, (λ, 1))
        else
            mult_eigs[end] = (λ, mult_eigs[end][2] + 1)
        end
    end
    
    return mult_eigs
end

@inline charpoly(A::AbstractMatrix) = berkowitz_vector(A) # Ordered from highest to lowest power
function berkowitz_vector(A::AbstractMatrix{T}) where {T}
    n = checksquare(A)
    @assert n >= 1

    M = T[1, -A[end, end]]

    @inbounds for i in 2:n
        # Loop from the back as that enables storing a single vector T,
        # rather than storing a matrix for each iteration before reducing.
        
        # i is the size of the A block to be partitioned into [a R; C Asub]

        Asub = @view(A[end - i + 2:end, end - i + 2:end])
        C = @view(A[end - i + 2:end, end - i + 1])      # Note: C is stored sequentially (it is a column)
        R = @view(A[end - i + 1, end - i + 2:end])      # Note R is not stored sequentially
        a = A[end - i + 1, end - i + 1]

        N = Vector{T}(undef, i + 1)
        N[1] = 1
        N[2] = -a

        K = -transpose(R)
        N[3] = K * C

        for j in 2:(i - 1)
            K *= Asub
            N[j + 2] = K * C
        end

        M = NonSquareLowerTriangularToeplitz(N, i) * M
    end

    return M
end

"""
Try to find symbolic roots of a characteristic polynomial with best effort. 
By Abel-Ruffini theorem, there is no general solution for polynomials of degree 5 or higher.
"""
function charpoly_roots(p::AbstractVector{T}) where {T <: IntOrRational}
    @assert length(p) >= 2 && isone(p[1])

    # p is ordered from highest to lowest power, hence if the last n coefficients are zero,
    # then we can divide through by x^n and we know that there will be n roots at 0.
    n = findfirst(i -> !iszero(p[i]), reverse(eachindex(p))) - 1
    eigs = zeros(T, n)

    p = p[1:end - n]
    
    if islinearorquadratic(p)
        return vcat(eigs, linearorquadraticroots(p))
    end

    # Only the leading and constant coefficients are non-zero. Take the nth root.
    if isbinomial(p)
        return vcat(eigs, binomialroots(p))
    end

    # Check if any of -5:5 are roots
    @inbounds for r in -5:5
        r = convert(T, r)
        while iszero(evalpoly(p, r))
            push!(eigs, r)
            p = divfactor(p, r)  # Divide out the factor (x - r)

            if islinearorquadratic(p)
                return vcat(eigs, linearorquadraticroots(p))
            end
        end
    end

    if iscubic(p)
        return vcat(eigs, cubicroots(p))
    end

    if isquartic(p)
        return vcat(eigs, quarticroots(p))
    end

    throw("Unable to find symbolic roots of polynomial")
end

islinear(p) = length(p) == 2
function linearroots(p)
    @assert islinear(p)

    a, b = p
    @assert isone(a)

    return [-b]
end

isquadratic(p) = length(p) == 3
function quadraticroots(p)
    @assert isquadratic(p)

    a, b, c = p
    @assert isone(a)

    if iszero(b)
        r = symbolic_sqrt(-c)
        return [r, -r]
    end

    d = b^2 - 4c
    if iszero(d)
        r = -symbolic_div(b, 2)
        return [r, r]
    else
        r1 = symbolic_div(-b + symbolic_sqrt(d), 2)
        r2 = symbolic_div(-b - symbolic_sqrt(d), 2)

        return [r1, r2]
    end
end

islinearorquadratic(p) = islinear(p) || isquadratic(p)
function linearorquadraticroots(p)
    if islinear(p)
        return linearroots(p)
    else
        return quadraticroots(p)
    end
end

isbinomial(p) = iszero(p[2:end - 1])
function binomialroots(p)
    @assert isbinomial(p)

    # TODO: Implement
    a, b = p[1], p[end]
    @assert isone(a)

    throw(NotImplementedError("Binomial polynomial roots is not implemented yet."))

    d = l - 1

    base = -b
    α = base^(1 // d)

    roots = [α * exp(k * 2π * im / d) for k in 0:(n - 1)]
    return roots
end

function evalpoly(p, r)
    coeff = reverse(p)
    f((i, c)) = c * r^(i - 1)
    return sum(f, enumerate(coeff))
end

function divfactor(p::AbstractVector{T}, r::T) where {T}
    n = length(p) - 1
    @assert n >= 1

    q_coeff = zeros(T, n)
    r_coeff = copy(p)

    @inbounds for i in 1:n
        q = r_coeff[i]
        q_coeff[i] = q

        r_coeff[i] -= q
        r_coeff[i + 1] -= -r * q
    end

    @assert iszero(r_coeff[end]) "r is not a factor of p"
    return q_coeff
end

iscubic(p) = length(p) == 4
function cubicroots(p)
    @assert iscubic(p)

    # Cardano's method
    a, b, c, d = p
    @assert isone(a)

    Q = symbolic_div(3c - b^2, 9)
    R = symbolic_div(9b * c - 27d - 2b^3, 54)

    D = Q^3 + R^2
    S = symbolic_cbrt(R + symbolic_sqrt(D))
    T = symbolic_cbrt(R - symbolic_sqrt(D))

    r1 = S + T - b // 3
    r2 = -symbolic_div(S + T, 2) - symbolic_div(b, 3) + (S - T) * symbolic_sqrt(3) * symbolic_div(1im, 2)
    r3 = -symbolic_div(S + T, 2) - symbolic_div(b, 3) - (S - T) * symbolic_sqrt(3) * symbolic_div(1im, 2)

    return [r1, r2, r3]
end

isquartic(p) = length(p) == 5
function quarticroots(p)
    @assert isquartic(p)

    a, b, c, d, e = p
    @assert isone(a)

    # Biquadratic equations
    if iszero(b) && iszero(d)
        r1, r2 = quadraticroots([1, c, e])
        return [symbolic_sqrt(r1), -symbolic_sqrt(r1), symbolic_sqrt(r2), -symbolic_sqrt(r2)]
    end

    # Ferrari's method
    q = [1, -c, b * d - 4e, 4c * e - b^2 * e - d^2]
    croots = cubicroots(q)
    y = findfirst(isreal, croots)

    p = [
        b + symbolic_sqrt(b^2 - 4c + 4y),
        b - symbolic_sqrt(b^2 - 4c + 4y)
    ]

    q = [
        y - symbolic_sqrt(y^2 - 4e),
        y + symbolic_sqrt(y^2 - 4e)
    ]

    r1 = symbolic_div(-p[1] + symbolic_sqrt(p[1]^2 - 8q[1]), 4)
    r2 = symbolic_div(-p[1] - symbolic_sqrt(p[1]^2 - 8q[1]), 4)
    r3 = symbolic_div(-p[2] + symbolic_sqrt(p[2]^2 - 8q[2]), 4)
    r4 = symbolic_div(-p[2] - symbolic_sqrt(p[2]^2 - 8q[2]), 4)

    return [r1, r2, r3, r4]
end

function symbolic_sqrt(r)
    if base == 0
        return 0
    end

    if r < 0
        return Symbolic.Pow(-r, 1//2) * 1im
    else
        return Symbolics.Pow(r, 1//2)
    end
end
symbolic_cbrt(r) = Symbolics.Pow(r, 1//3)
symbolic_div(num, den) = Symbolics.div(num, den)