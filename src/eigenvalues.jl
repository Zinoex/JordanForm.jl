function radical_eigenvalues(A::AbstractMatrix{T}) where {T <: IntOrRational}
    n = checksquare(A)
    A = convert.(Rational{eltype(A)}, A)
    
    p = charpoly(A)
end

function charpoly(A::AbstractMatrix{<:Rational})
    n = checksquare(A)

    p = berkowitz_vector(A, n)

    return p
end


function berkowitz_vector(A::AbstractMatrix{T}, n) where {T}
    M = ones(T, 1)

    for i in 1:(n - 1)  
        # Loop form the back as that enables storing a single vector T,
        # rather than storing a matrix for each iterations before reducing.

        if i == 1
            M = T[1; -A[end, end]] * M
        elseif i == 2
            C = A[end, end - 1]
            R = A[end - 1, end]
            a = A[end - 1, end - 1]

            M = NonSquareLowerTriangularToeplitz(T[1, -a,  -R * C], i + 1) * M
        else
            Asub = @view(A[end - i:end, end - i:end])
            C = @view(A[end - i:end, end - i - 1])      # Note: C is stored sequentially (it is a column)
            R = @view(A[end - i - 1, end - i:end])      # Note R is not stored sequentially
            a = A[end - i - 1, end - i - 1]

            N = Vector{T}(undef, i + 1)
            N[1] = 1
            N[2] = -a

            K = -R
            N[3] = K * C     # Multiplying by C is faster than multiplying by R due to the storage. Hence we multiply by C in the iteration.

            for j in 2:(i - 1)
                K *= Asub
                N[j + 2] = K * C
            end

            M = NonSquareLowerTriangularToeplitz(N) * M
        end
    end

    return M
end