using IterTools, LinearAlgebra, Requires, SparseArrays

using MultivariatePolynomials: AbstractVariable, AbstractMonomialLike, exponents, variables, powers,
                               monomials, coefficient
using ReachabilityBase.Arrays: logarithmic_norm

function __init__()
    @require LazySets = "b4f0291d-fe17-52bc-9479-3d1a343d9043" begin
        using .LazySets: Hyperrectangle, dim, low, high
    end
    @require IntervalArithmetic = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253" begin
        using .IntervalArithmetic: interval
    end
end
