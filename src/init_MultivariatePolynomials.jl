export quadratic_matrix_form

eval(quote
    using .MultivariatePolynomials: AbstractVariable, AbstractMonomialLike,
                                    exponents, variables, powers, monomials, coefficient
end)

eval(load_kron_multivariate())
eval(load_quadratic_matrix_form())
