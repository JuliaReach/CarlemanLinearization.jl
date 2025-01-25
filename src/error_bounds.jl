# =========================================================================================================================
# Theoretical error bounds for Carleman linearization
#
# References:
#
# - [ForetsP17](@citet)
#
# - [LiuKKLTC21](@citet)
#
# =========================================================================================================================

# --- Error bounds using a priori estimate from [ForetsP17](@citet) ---
# These bounds use the supremum norm (p = Inf).

# See [ForetsP17; Theorem 4.2](@citet); this bound is based on an a priori
# estimate of the norm of the exact solution x(t).
function error_bound_apriori(α, F₁, F₂; N)
    nF₂ = opnorm(F₂, Inf)
    μF₁ = logarithmic_norm(F₁, Inf)

    β = α * nF₂ / μF₁
    ε = t -> α * β^N * (exp(μF₁ * t) - 1)^N
    return ε
end

# See [ForetsP17; Theorem 4.2](@citet).
function convergence_radius_apriori(α, F₁, F₂)
    nF₂ = opnorm(F₂, Inf)
    μF₁ = logarithmic_norm(F₁, Inf)

    if μF₁ < 0
        return Inf
    end
    β = α * nF₂ / μF₁
    T = (1 / μF₁) * log(1 + 1 / β)
    return T
end

# --- Error bounds using power series method from [ForetsP17](@citet) ---

# See [ForetsP17; Theorem 4.3](@citet), which uses the power series method.
function error_bound_pseries(x₀, F₁, F₂; N)
    nx₀ = norm(x₀, Inf)
    nF₁ = opnorm(F₁, Inf)
    nF₂ = opnorm(F₂, Inf)
    β₀ = nx₀ * nF₂ / nF₁

    ε = t -> nx₀ * exp(nF₁ * t) / (1 - β₀ * (exp(nF₁ * t) - 1)) * (β₀ * (exp(nF₁ * t) - 1))^N
    return ε
end

# See [ForetsP17; Theorem 4.3](@citet).
function convergence_radius_pseries(x₀, F₁, F₂)
    nx₀ = norm(x₀, Inf)
    nF₁ = opnorm(F₁, Inf)
    nF₂ = opnorm(F₂, Inf)
    β₀ = nx₀ * nF₂ / nF₁

    T = (1 / nF₁) * log(1 + 1 / β₀)
    return T
end

# --- Error bounds using spectral abscissa from [LiuKKLTC21](@citet) ---
# These bounds use the spectral norm (p = 2).

# compute eigenvalues and sort them by increasing real part
function _error_bound_specabs_Re_λ₁(F₁)
    λ = eigvals(F₁; sortby=real)
    return real(last(λ))
end

# See [LiuKKLTC21; Corollary 1](@citet).
function error_bound_specabs(x₀, F₁, F₂; N, check=true)
    Re_λ₁ = _error_bound_specabs_Re_λ₁(F₁)
    if check && Re_λ₁ >= 0
        throw(ArgumentError("expected Re(λ₁) < 0, got $Re_λ₁"))
    end

    # Equation (2.2)
    nx₀ = norm(x₀, 2)
    nF₂ = opnorm(F₂, 2)
    R = nx₀ * nF₂ / -Re_λ₁
    if check && R >= 1
        throw(ArgumentError("expected R < 1, got R = $R; try scaling the ODE"))
    end

    # Equation (4.29)
    ε = t -> nx₀ * R^N * (1 - exp(Re_λ₁ * t))^N
    return ε
end

# See [LiuKKLTC21; Corollary 1](@citet).
function convergence_radius_specabs(x₀, F₁, F₂)
    Re_λ₁ = _error_bound_specabs_Re_λ₁(F₁)

    if Re_λ₁ < 0
        T = Inf
    elseif iszero(Re_λ₁)
        nx₀ = norm(x₀, 2)
        nF₂ = opnorm(F₂, 2)
        β = nx₀ * nF₂
        T = 1 / β
    else
        throw(ArgumentError("expected spectral abscissa to be negative or zero, got $Re_λ₁"))
    end
    return T
end
