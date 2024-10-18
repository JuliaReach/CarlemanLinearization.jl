using Test, CarlemanLinearization
import Aqua

@testset "Aqua tests" begin
    Aqua.test_all(CarlemanLinearization)
end
