using LinearAlgebra: /, Diagonal, eigvals, kron, norm, opnorm
using SparseArrays: findnz, sparse, spzeros
using Requires: @require

using MultivariatePolynomials: AbstractVariable, monomials, coefficient
using ReachabilityBase.Arrays: logarithmic_norm
using ReachabilityBase.Require: require

function __init__()
    @require LazySets = "b4f0291d-fe17-52bc-9479-3d1a343d9043" begin
        using .LazySets: Hyperrectangle, dim, low, high
    end
    @require IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253" begin
        using .IntervalArithmetic: interval, inf, sup
    end
end
