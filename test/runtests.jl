using Test, CarlemanLinearization
using DynamicPolynomials, MultivariatePolynomials

using LazySets: Hyperrectangle, low, high
using CarlemanLinearization: generate_monomials, build_matrix, kron_pow,
                             quadratic_matrix_form, lift_vector

@testset "Kronecker power (symbolic)" begin
    @polyvar x[1:2]
    y = kron_pow(x, 2)
    @test findall(==(x[1] * x[2]), y) == [2, 3]
    @test findall(==(x[1]^2), y) == [1]
    @test findall(==(x[2]^2), y) == [4]
    @test findall(==(x[2]^3), y) == Int[]
end

@testset "Conversion from polynomial to matrix representation" begin
    vars = @polyvar x y
    dx = 3x + y^2
    dy = x - y - 2.2 * x * y + x^2
    f = [dx, dy]
    F1, F2 = quadratic_matrix_form(f, vars)
    @test F1 == [3.0 0; 1 -1]
    @test F2 == [0.0 0 0 1; 1 -2.2 0 0]
end

@testset "Generating all commutative monomials" begin
    for (n, D) in [(3, 0), (3, 1), (3, 2), (5, 3), (5, 10)]
        monomials = generate_monomials(n, D)
        # No repetitions
        @test length(monomials) == length(Set(monomials))
        # Correct size
        @test length(monomials) == binomial(n + D, n)
        # Nonnegativity
        @test all([sum(m) == sum(map(abs, m)) for m in monomials])
        # ordering
        @test all([sum(monomials[i - 1]) <= sum(monomials[i]) for i in 2:length(monomials)])
        # not exceeding D
        @test sum(monomials[end]) == D
    end
end

@testset "Compressed matrices" begin
    F1 = [3.0 0; 1 -1]
    F2 = [0.0 0 0 1; 1 -2.2 0 0]
    result = build_matrix(F1, F2, 1; compress=true)
    @test result == F1

    result = build_matrix(F1, F2, 2; compress=true)
    @test result == [3.0 0.0 0.0 0.0 1.0;
                     1.0 -1.0 1.0 -2.2 0;
                     0 0 6.0 0 0;
                     0 0 1.0 -1.0 0;
                     0 0 0 2.0 -2.0]
end

@testset "Lifting vectors" begin
    x0 = Hyperrectangle(; low=[0, 1, -1], high=[1, 2, 3])
    N = 2
    lifted = lift_vector(x0, N)
    sides = Set(zip(low(lifted), high(lifted)))
    @test sides ==
          Set{Tuple{Rational,Rational}}([(0, 1), (1, 2), (-1, 3), (0, 1), (1, 4), (0, 9), (0, 2),
                                         (-1, 3), (-2, 6)])
end
