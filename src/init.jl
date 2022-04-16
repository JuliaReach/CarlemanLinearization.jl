using IterTools, LinearAlgebra, SparseArrays

using MultivariatePolynomials: AbstractVariable, AbstractMonomialLike,
      exponents, variables, powers, monomials, coefficient

using LazySets: Hyperrectangle, dim, low, high
using IntervalArithmetic: interval
