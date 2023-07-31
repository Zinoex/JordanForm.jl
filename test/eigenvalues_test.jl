@testset "Cubic diagonal" begin
    λ_1 = 1
    λ_2 = 2
    λ_3 = 3
    A = [λ_1 0 0; 0 λ_2 0; 0 0 λ_3]
    eigs = radical_eigvals(A)
    @test all(eigs .=== [λ_1, λ_2, λ_3])
    @test all(algebraic_multiplicity(eigs) .=== [(λ_1, 1), (λ_2, 1), (λ_3, 1)])
end

@testset "Quartic diagonal" begin
    λ_1 = 1
    λ_2 = 2
    λ_3 = 3
    λ_4 = 4
    A = [λ_1 0 0 0; 0 λ_2 0 0; 0 0 λ_3 0; 0 0 0 λ_4]
    eigs = radical_eigvals(A)
    @test all(eigs .=== [λ_1, λ_2, λ_3, λ_4])
    @test all(algebraic_multiplicity(eigs) .=== [(λ_1, 1), (λ_2, 1), (λ_3, 1), (λ_4, 1)])
end

@testset "Cubic not diagonal" begin
    λ_1 = 1
    λ_2 = 2
    λ_3 = 3
    A = [1 1 0; 0 2 1; 0 0 3]
    eigs = radical_eigvals(A)
    @test all(eigs .=== [λ_1, λ_2, λ_3])
    @test all(algebraic_multiplicity(eigs) .=== [(λ_1, 1), (λ_2, 1), (λ_3, 1)])
end
