export quadratic_matrix_form

eval(quote
    using .DynamicPolynomials: Monomial, PolyVar
end)

eval(load_quadratic_matrix_form())
