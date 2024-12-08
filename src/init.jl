using LinearAlgebra, SparseArrays

using MultivariatePolynomials: AbstractVariable, AbstractMonomialLike, exponents, variables, powers,
                               monomials, coefficient
using ReachabilityBase.Arrays: logarithmic_norm

# functionality available through optional dependencies
using PackageExtensionCompat
function __init__()
    @require_extensions
end
function lift_vector end
