using Test, CarlemanLinearization
import Aqua

import Pkg
@static if VERSION >= v"1.6"  # TODO make explicit test requirement
    Pkg.add("ExplicitImports")
    import ExplicitImports

    @testset "ExplicitImports tests" begin
        ignores = (:generate_monomials,)  # false positive due to package extensions
        @test isnothing(ExplicitImports.check_all_explicit_imports_are_public(CarlemanLinearization;
                                                                              ignore=ignores))
        @test isnothing(ExplicitImports.check_all_explicit_imports_via_owners(CarlemanLinearization))
        @test isnothing(ExplicitImports.check_all_qualified_accesses_are_public(CarlemanLinearization))
        @test isnothing(ExplicitImports.check_all_qualified_accesses_via_owners(CarlemanLinearization))
        @test isnothing(ExplicitImports.check_no_implicit_imports(CarlemanLinearization))
        @test isnothing(ExplicitImports.check_no_self_qualified_accesses(CarlemanLinearization))
        @test isnothing(ExplicitImports.check_no_stale_explicit_imports(CarlemanLinearization))
    end
end

@static if VERSION >= v"1.10"
    # JET v0.9.0 (earliest supported version) requires Julia v1.10
    Pkg.add("JET")
    import JET

    @testset "JET tests" begin
        # false positives for Base functionality
        JET.test_package(CarlemanLinearization; target_modules=(CarlemanLinearization,))
    end
end

@testset "Aqua tests" begin
    Aqua.test_all(CarlemanLinearization)
end
