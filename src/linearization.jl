# ======================================
# Construction of the linearized system
# ======================================

"""
    build_matrix(F₁, F₂, N)

Compute the Carleman linearization matrix associated to the quadratic
system ``x' = F₁x + F₂(x⊗x)``, truncated at order ``N``.

### Input

- `F₁` -- sparse matrix of size ``n × n``
- `F₂` -- sparse matrix of size ``n × n^2``
- `N`  -- integer representing the truncation order, should be at least two

### Output

Sparse matrix `A`.

### Algorithm

See references [1] and [2] of CARLIN.md.
"""
function build_matrix(F₁, F₂, N)
    if N < 2
        throw(ArgumentError("expected the truncation order to be at least 2, got N=$N"))
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

function _build_matrix(F₁, F₂, N)
    @assert N >= 3 "expected N to be at least 3, got N=$N"
    n = size(F₁, 1)

    out = Vector{typeof(F₁)}()
    for j in 1:N
        if j == N
            F12_j = kron_sum(F₁, j)
            Zleft_j = spzeros(n^j, sum(n^i for i in 1:j-1))
            push!(out, hcat(Zleft_j, F12_j))
        else
            F12_j = hcat(kron_sum(F₁, j), kron_sum(F₂, j))
            if j == 1
                Zright_j = spzeros(n^j, sum(n^i for i in j+2:N))
                push!(out, hcat(F12_j, Zright_j))
            elseif j == N-1
                Zleft_j = spzeros(n^j, sum(n^i for i in 1:j-1))
                push!(out, hcat(Zleft_j, F12_j))
            else # general case
                Zleft_j = spzeros(n^j, sum(n^i for i in 1:j-1))
                Zright_j = spzeros(n^j, sum(n^i for i in j+2:N))
                push!(out, hcat(Zleft_j, F12_j, Zright_j))
            end
        end
    end
    return reduce(vcat, out)
end
