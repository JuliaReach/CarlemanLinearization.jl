# ======================
# Hard dependencies
# ======================

using LinearAlgebra, SparseArrays, Requires, IterTools

using LazySets, IntervalArithmetic

using ReachabilityAnalysis
using LazySets.Arrays

# functions for Kronecker powers and sums
#using ReachabilityAnalysis: kron_pow, kron_pow_stack
using ReachabilityAnalysis: ReachSolution, LGG09


# ======================
# Optional dependencies
# ======================

using Requires

function __init__()
    # tools for symbolic algebra
    @require MultivariatePolynomials = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3" include("init_MultivariatePolynomials.jl")
end
