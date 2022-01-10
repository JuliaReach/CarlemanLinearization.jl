# ======================
# Hard dependencies
# ======================

using LinearAlgebra, SparseArrays, Requires

# ======================
# Optional dependencies
# ======================

using Requires

function __init__()
    # tools for symbolic algebra
    @require MultivariatePolynomials = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3" include("init_MultivariatePolynomials.jl")

    # sparse dynamic representation of multivariate polynomials
    @require DynamicPolynomials = "7c1d4256-1411-5781-91ec-d7bc3513ac07" include("init_DynamicPolynomials.jl")
end
