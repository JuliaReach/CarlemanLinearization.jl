using LinearAlgebra, SparseArrays

using MultivariatePolynomials: AbstractVariable, AbstractMonomialLike, exponents, variables, powers,
                               monomials, coefficient
using ReachabilityBase.Arrays: logarithmic_norm
using ReachabilityBase.Require

# functionality available through optional dependencies
using PackageExtensionCompat
function __init__()
    @require_extensions
end
