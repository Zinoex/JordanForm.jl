function radical_eigenvalues(A::AbstractMatrix{T}) where {T <: IntOrRational}
    n = LinearAlgebra.checksquare(A)
    
    p = charpoly(A)
end

function charpoly(A::AbstractMatrix)
    p = berkowitz_vector(A)

    return p
end


function berkowitz_vector(A::AbstractMatrix{T}) where {T}
    n = LinearAlgebra.checksquare(A)
    @assert n >= 1

    M = T[1, -A[end, end]]

    for i in 2:n
        # Loop form the back as that enables storing a single vector T,
        # rather than storing a matrix for each iterations before reducing.
        
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