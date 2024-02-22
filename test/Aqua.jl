using Test, CarlemanLinearization
import Aqua

@testset "Aqua tests" begin
    Aqua.test_all(CarlemanLinearization; ambiguities=false)

    # do not warn about ambiguities in dependencies
    Aqua.test_ambiguities([CarlemanLinearization, Core])
end
