using Base: /, kron
using LinearAlgebra: Diagonal, eigvals, norm, opnorm
using SparseArrays: findnz, sparse, spzeros

using MultivariatePolynomials: AbstractVariable, monomials, coefficient
using ReachabilityBase.Arrays: logarithmic_norm
using ReachabilityBase.Require: require

@static if !isdefined(Base, :get_extension)
    using Requires: @require
    # COV_EXCL_START
    function __init__()
        @require LazySets = "b4f0291d-fe17-52bc-9479-3d1a343d9043" begin
            @require IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253" begin
                include("../ext/LazySetsExt.jl")
            end
        end
    end
    # COV_EXCL_STOP
end
