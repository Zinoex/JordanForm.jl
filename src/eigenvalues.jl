function radical_eigvals(A::AbstractMatrix{T}) where {T <: IntOrRational}
    LinearAlgebra.checksquare(A)
    
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
end

@inline charpoly(A::AbstractMatrix) = berkowitz_vector(A) # Ordered from highest to lowest power
function berkowitz_vector(A::AbstractMatrix{T}) where {T}
    n = LinearAlgebra.checksquare(A)
    @assert n >= 1

    M = T[1, -A[end, end]]

    for i in 2:n
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

function charpoly_roots(p::AbstractVector{T}) where {T <: IntOrRational}
    @assert length(p) >= 2 && isone(p[1])

    # p is ordered from highest to lowest power, hence if the last n coefficients are zero,
    # then we know that there will be n roots at 0, and we can divide through by x^n.
    n = findfirst(i -> !iszero(p[i]), n:1) - 1
    eigs = zeros(T, n)

    p = p[1:end - n]
    l = length(p)
    
    if l == 2  # It is a linear polynomial
        push!(eigs, -p[2])
        return eigs
    elseif l == 3  # It is a quadratic polynomial
        a, b, c = p

        if iszero(b)
            r = Symbolics.wrap(Symbolics.Term(sqrt, -c // a))
            push!(eigs, r)
            push!(eigs, -r)
            return eigs
        end

        d = b^2 - 4a*c
        if iszero(d)
            push!(eigs, -b // (2a))
            push!(eigs, -b // (2a))
        else
            push!(eigs, (-b + sqrt(d)) / (2a))
            push!(eigs, (-b - sqrt(d)) / (2a))
        end

        return eigs
    end

    # Only the leading and constant coefficients are non-zero. Take the nth root.
    if iszero(p[2:end - 1])
        a, b = p[1], p[end]
        d = l - 1

        base = -b
        α = base^(1 // d)

        roots = [α * exp(k * 2π * im / d) for k in 0:(n - 1)]
        return vcat(eigs, roots)
    end

    if l == 4
        # TODO: Handle the quartic polynomial case
    end

    # TODO: Better, more meaningful error message
    throw(NotImplementedError("The charpoly_roots function has not been implemented for polynomials of degree $l"))
end

function symbolic_radical(base, n)
    if base == 0
        return 0
    end

    if base < 0
        imagimary = 1im
        base = -base
    else
        return base^(1 / n)
    end
end