using Test, CarlemanLinearization
using DynamicPolynomials, MultivariatePolynomials, LinearAlgebra

using LazySets: Hyperrectangle, low, high
using CarlemanLinearization: generate_monomials, _build_matrix_N

@testset "Kronecker power of identity matrix" begin
    @test kron_id(2, 3) == Diagonal(ones(8))
end

@testset "Kronecker power (symbolic)" begin
    @polyvar x[1:2]
    y = kron_pow(x, 2)
    @test findall(==(x[1] * x[2]), y) == [2, 3]
    @test findall(==(x[1]^2), y) == [1]
    @test findall(==(x[2]^2), y) == [4]
    @test findall(==(x[2]^3), y) == Int[]
end

@testset "Kronecker sandwich I^{⊗ l} ⊗ F ⊗ I^{⊗ r}" begin
    F = [0 1; -1 0.0]
    I2 = I(2)

    Q = kron_sandwich(F, 2, 2, 2)
    @test size(Q) == (32, 32)
    @test reduce(kron, [I2, I2, F, I2, I2]) == Q

    Q = kron_sandwich(F, 2, 0, 2)
    @test size(Q) == (8, 8)
    @test reduce(kron, [F, I2, I2]) == Q

    Q = kron_sandwich(F, 2, 2, 0)
    @test size(Q) == (8, 8)
    @test reduce(kron, [I2, I2, F]) == Q
end

@testset "Kronecker sum" begin
    F = [0 1; -1 0.0]
    I2 = I(2)
    @test_throws ArgumentError kron_sum(F, 0)
    @test kron_sum(F, 1) == F
    @test kron_sum(F, 2) == kron(F, I2) + kron(I2, F)
    @test kron_sum(F, 3) ==
          sum(reduce(kron, X) for X in [[F, I2, I2], [I2, F, I2], [I2, I2, F]])
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

    @test_throws ArgumentError build_matrix(F1, F2, 0; compress=true)

    result = build_matrix(F1, F2, 2; compress=true)
    @test result == [3.0 0.0 0.0 0.0 1.0;
                     1.0 -1.0 1.0 -2.2 0;
                     0 0 6.0 0 0;
                     0 0 1.0 -1.0 0;
                     0 0 0 2.0 -2.0]

    result = build_matrix(F1, F2, 2; compress=false)
    @test result == [kron_sum(F1, 1) kron_sum(F2, 1);
                     zeros(4, 2) kron_sum(F1, 2)]

    result = build_matrix(F1, F2, 3; compress=false)
    @test result == [kron_sum(F1, 1) kron_sum(F2, 1) zeros(2, 8);
                     zeros(4, 2) kron_sum(F1, 2) kron_sum(F2, 2);
                     zeros(8, 6) kron_sum(F1, 3)]
    @test result == _build_matrix_N(F1, F2, 3)

    result = build_matrix(F1, F2, 4; compress=false)
    @test result == [kron_sum(F1, 1) kron_sum(F2, 1) zeros(2, 24);
                     zeros(4, 2) kron_sum(F1, 2) kron_sum(F2, 2) zeros(4, 16);
                     zeros(8, 6) kron_sum(F1, 3) kron_sum(F2, 3);
                     zeros(16, 14) kron_sum(F1, 4)]
    @test result == _build_matrix_N(F1, F2, 4)

    result = build_matrix(F1, F2, 5; compress=false)
    @test size(result) == (62, 62)
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

@testset "Error bounds" begin
    F1 = [0.0 1; -1 -2]
    F2 = [0 0 0 1; 1 -2.2 0 0]
    x0 = 0.1

    # error bounds using power-series method

    T = convergence_radius_pseries(x0, F1, F2)
    @test T ≈ 0.7797 atol = 1e-4

    e2 = error_bound_pseries(x0, F1, F2; N=2)
    e3 = error_bound_pseries(x0, F1, F2; N=3)
    for t in 0.01:0.01:T
        @test e3(t) < e2(t)
    end

    # error bounds using spectral abscissa

    # convergence everywhere (T = Inf)
    T = convergence_radius_specabs(x0, F1, F2)
    @test T == Inf

    e2 = error_bound_specabs(x0, F1, F2; N=2)
    e3 = error_bound_specabs(x0, F1, F2; N=3)
    for t in 0.01:0.01:1
        @test e3(t) < e2(t)
    end

    # R >= 1
    @test_throws ArgumentError error_bound_specabs(5 * x0, F1, F2; N=1)

    # Re_λ₁ = 0
    F1_zero = [0.0 0; 0 -1]
    T = convergence_radius_specabs(x0, F1_zero, F2)

    @test_throws ArgumentError error_bound_specabs(x0, F1_zero, F2; N=2)

    # Re_λ₁ > 0
    F1_pos = [0.0 0; 0 1]
    @test_throws ArgumentError convergence_radius_specabs(x0, F1_pos, F2)
    @test_throws ArgumentError convergence_radius_specabs(x0, F1_pos, F2; check=false)
end
