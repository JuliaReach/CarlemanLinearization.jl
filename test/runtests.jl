using Test, CarlemanLinearization
using DynamicPolynomials

@testset "Kronecker power (symbolic)" begin
    @polyvar x[1:2]
    y = kron_pow(x, 2)
    @test findall(y, x[1]*x[2]) == [2, 3]
    @test findall(y, x[1]^2) == [1]
    @test findall(y, x[2]^2) == [4]
    @test findall(y, x[2]^3) == Int[]
end
