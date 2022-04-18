using IterTools, LinearAlgebra, Requires, SparseArrays

using MultivariatePolynomials: AbstractVariable, AbstractMonomialLike,
      exponents, variables, powers, monomials, coefficient

function __init__()
    @require LazySets = "b4f0291d-fe17-52bc-9479-3d1a343d9043" begin
        using .LazySets: Hyperrectangle, dim, low, high
        using .LazySets.IntervalArithmetic: interval
    end
end
