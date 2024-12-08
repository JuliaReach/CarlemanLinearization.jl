# ======================================
# Construction of the linearized system
# ======================================

"""
    build_matrix(F₁, F₂, N; compress=false)

Compute the Carleman linearization matrix associated to the quadratic
system ``x' = F₁x + F₂(x⊗x)``, truncated at order ``N``.

### Input

- `F₁` -- sparse matrix of size ``n × n``
- `F₂` -- sparse matrix of size ``n × n^2``
- `N`  -- integer representing the truncation order, should be at least two
- `compress` -- if `true` produces the matrix corresponding to the commutative monomials only

### Output

Sparse matrix `A`.

### Algorithm

See references [1] and [2] of CARLIN.md.
"""
function build_matrix(F₁, F₂, N; compress=false)
    F₁ = sparse(F₁)
    F₂ = sparse(F₂)
    if N > 1 && compress
        return build_matrix_compressed(F₁, F₂, N)
    end

    # No compression
    if N < 1
        throw(ArgumentError("expected the truncation order to be at least 1, got N=$N"))
    elseif N == 1
        _build_matrix_N1(F₁, F₂)
    elseif N == 2
        _build_matrix_N2(F₁, F₂)
    elseif N == 3
        _build_matrix_N3(F₁, F₂)
    elseif N == 4
        _build_matrix_N4(F₁, F₂)
    else
        _build_matrix_N(F₁, F₂, N) # general case
    end
end

"""
    build_matrix_compressed(F₁, F₂, N)

Compute the compressed Carleman linearization matrix associated to the quadratic
system ``x' = F₁x + F₂(x⊗x)``, truncated at order ``N``.

Input & Output are the same as for `build_matrix`.
"""
function build_matrix_compressed(F₁, F₂, N)
    n = size(F₁)[1]
    monoms = generate_monomials(n, N)
    # skip the first monomial, which is always the constant 1
    nonfirst_monoms = Iterators.peel(monoms)[2]
    monom_to_ind = Dict(m => i for (i, m) in enumerate(nonfirst_monoms))
    result = spzeros(length(monoms) - 1, length(monoms) - 1)
    # linear/quadratic monomials on the right-hand side
    linear_rhs, quadratic_rhs = findnz(F₁), findnz(F₂)
    for (ind, m) in enumerate(nonfirst_monoms)
        # for a given monomial m of degree d, we compute the degree-d part of
        # its derivative m' using the linear part, F₁, of the ode system
        for (i, j, c) in zip(linear_rhs...)
            if m[i] > 0
                deriv = m
                if i != j
                    deriv = m .+ Tuple((k == i) ? -1 : ((k == j) ? 1 : 0) for k in 1:n)
                end
                result[ind, monom_to_ind[deriv]] = m[i] * c
            end
        end

        # for a given monomial m of degree d, we compute the degree-(d + 1) part
        # of its derivative m' using the linear part, F₂, of the ode system
        if sum(m) < N
            for (i, j, c) in zip(quadratic_rhs...)
                # extracting the indices of the variables corresponding to the
                # degree-2 monomial corresponding to the j-th column of F₂
                j0 = ((j - 1) % n) + 1
                j1 = ((j - 1) ÷ n) + 1
                if m[i] > 0
                    deriv = [m...]
                    deriv[i] -= 1
                    deriv[j0] += 1
                    deriv[j1] += 1
                    deriv = (deriv...,)
                    result[ind, monom_to_ind[deriv]] += m[i] * c
                end
            end
        end
    end
    return result
end

"""
    generate_monomials(n, N)

Return a list of `n`-tuples of nonegative integers with the sum at most `N`,
ordered by the total degree (no other guarantees on the ordering).
"""
function generate_monomials(n, N)
    if n == 1
        return [(i,) for i in 0:N]
    end
    result = []
    prev = generate_monomials(n - 1, N)
    for ord in 0:N
        pind = 1
        while pind <= length(prev) && sum(prev[pind]) <= ord
            push!(result, (ord - sum(prev[pind]), prev[pind]...))
            pind += 1
        end
    end
    return result
end

function _build_matrix_N1(F₁, _)
    return F₁
end

function _build_matrix_N2(F₁, F₂)
    n = size(F₁, 1)
    a = hcat(kron_sum(F₁, 1), kron_sum(F₂, 1))
    b = hcat(spzeros(n^2, n), kron_sum(F₁, 2))
    return vcat(a, b)
end

function _build_matrix_N3(F₁, F₂)
    n = size(F₁, 1)
    a = hcat(kron_sum(F₁, 1), kron_sum(F₂, 1), spzeros(n, n^3))
    b = hcat(spzeros(n^2, n), kron_sum(F₁, 2), kron_sum(F₂, 2))
    c = hcat(spzeros(n^3, n), spzeros(n^3, n^2), kron_sum(F₁, 3))
    return vcat(a, b, c)
end

function _build_matrix_N4(F₁, F₂)
    n = size(F₁, 1)
    a = hcat(kron_sum(F₁, 1), kron_sum(F₂, 1), spzeros(n, n^3), spzeros(n, n^4))
    b = hcat(spzeros(n^2, n), kron_sum(F₁, 2), kron_sum(F₂, 2), spzeros(n^2, n^4))
    c = hcat(spzeros(n^3, n), spzeros(n^3, n^2), kron_sum(F₁, 3), kron_sum(F₂, 3))
    d = hcat(spzeros(n^4, n), spzeros(n^4, n^2), spzeros(n^4, n^3), kron_sum(F₁, 4))
    return vcat(a, b, c, d)
end

function _build_matrix_N(F₁, F₂, N)
    @assert N >= 3 "expected N to be at least 3, got N = $N"
    n = size(F₁, 1)

    out = Vector{typeof(F₁)}(undef, N)
    for j in 1:N
        if j == N
            F12_j = kron_sum(F₁, j)
            Zleft_j = spzeros(n^j, sum(n^i for i in 1:(j - 1)))
            out[j] = hcat(Zleft_j, F12_j)
        else
            F12_j = hcat(kron_sum(F₁, j), kron_sum(F₂, j))
            if j == 1
                Zright_j = spzeros(n^j, sum(n^i for i in (j + 2):N))
                out[j] = hcat(F12_j, Zright_j)
            elseif j == N - 1
                Zleft_j = spzeros(n^j, sum(n^i for i in 1:(j - 1)))
                out[j] = hcat(Zleft_j, F12_j)
            else # general case
                Zleft_j = spzeros(n^j, sum(n^i for i in 1:(j - 1)))
                Zright_j = spzeros(n^j, sum(n^i for i in (j + 2):N))
                out[j] = hcat(Zleft_j, F12_j, Zright_j)
            end
        end
    end
    return reduce(vcat, out)
end
