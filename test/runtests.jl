using Test, CarlemanLinearization
using DynamicPolynomials

@testset "Kronecker power (symbolic)" begin
    @polyvar x[1:2]
    y = kron_pow(x, 2)
    @test findall(==(x[1]*x[2]), y) == [2, 3]
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
