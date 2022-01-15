function load_quadratic_matrix_form()
return quote

"""
    quadratic_matrix_form(f, vars)

Extract the linear and quadratic matrices of a system of polynomial equations.

### Input

- `f`    -- vector of polynomials
- `vars` -- vector of variables

### Output

Matrices `F1` and `F2` of size `n × n` and `n × n^2` respectively,
such that `F1` contains the linear coefficients of the polynomial vector
field `f` and `F2` the quadratic coefficients, which indexing that respects
the Kronecker power `vars ⊗ vars`.

### Examples

TO-DO
"""
function quadratic_matrix_form(f, vars)
    n = length(vars)
    @assert length(f) == n

    # get the linear part
    F1 = zeros(n, n)
    update!(F1, vars, f)

    # get the quadratic part
    F2 = zeros(n, n^2)
    update!(F2, kron_pow(vars, 2), f)

    return F1, F2
end

function update!(Fk, xk, f)
    n = size(Fk, 1)
    for i in 1:n
        monom = monomials(f[i])
        S = Set()
        for (j, x) in enumerate(xk)
            # skip monomials already seen
            x ∈ S ? continue : push!(S, x)
            idx = findfirst(==(x), monom)
            if !isnothing(idx)
                Fk[i, j] = coefficient(f[i], x)
            end
        end
    end
    return Fk
end

end end  # quote / load_quadratic_matrix_form()
