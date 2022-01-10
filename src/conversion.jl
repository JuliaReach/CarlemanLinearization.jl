function load_quadratic_matrix_form()
return quote

function _get_vars_indices(eqs)
    V = typeof(eqs[1][1])
    x2i = Dict{V, Int}()
    vars = Vector{V}(undef, length(eqs))
    for i in eachindex(eqs)
        x = eqs[i][1]
        x in keys(x2i) && error("duplicate equation for variable `$x`")
        vars[i] = x
        x2i[x] = i
    end
    return vars, x2i
end

function quadratic_matrix_form(eqs::AbstractVector)
    vars, x2i = _get_vars_indices(eqs)
    extvars = kron_pow(vars, 2)
    # add all variables to each monomial (needed for compatibility later)
    n = length(vars)
    for (i, ti) in enumerate(extvars)
        z = Vector{Int}(undef, n)
        k = 1
        for (j, vj) in enumerate(vars)
            if vj in ti.vars
                z[j] = ti.z[k]  # assumes the same variable order
                k += 1
            else
                z[j] = 0
            end
        end
        extvars[i] = Monomial(vars, z)
    end
    # replace duplicates by a fresh monomial (needed for compatibility later)
    _replace_duplicates!(extvars)

    n = length(eqs)
    F1 = zeros(n, n)
    F2 = zeros(n, n^2)

    @inbounds for (i, (x, dx)) in enumerate(eqs)
        quadratic_matrix_form!(F1, F2, x, dx, x2i=x2i, vars=vars, extvars=extvars)
    end
    return F1, F2
end

function _replace_duplicates!(x)
    fresh = Monomial(PolyVar{true}("___dummy___"))  # assumed to not exist
    seen = Set()
    @inbounds for (i, vi) in enumerate(x)
        if vi in seen
            x[i] = fresh
        else
            push!(seen, vi)
        end
    end
end

function _onehot(i, n)
    v = zeros(Int, n)
    v[i] = 1
    return v
end

function quadratic_matrix_form!(F1, F2, x, dx; x2i, vars, extvars)
    row = x2i[x]

    # linear monomials (containing all variables for compatibility later)
    n = length(vars)
    lins = [Monomial(vars, _onehot(i, n)) for i in eachindex(vars)]

    # find summands in extvars
    @inbounds for summand in dx
        polyn = DynamicPolynomials.polynomial(summand)
        monom = DynamicPolynomials.monomial(summand.x)
        coeff = Float64(polyn.a[1])
        j = findfirst(lins, monom)
        if !isnothing(j)
            # linear term
            F1[row, j] = coeff
            continue
        end
        j = findall(extvars, monom)
        if !isempty(j)
            # quadratic term
            for jj in j
                F2[row, jj] = coeff
            end
            continue
        end
    end
end

end end  # quote / load_quadratic_matrix_form()
