# =========================================================================================================================
# Theoretical error bounds for Carleman linearization
#
# References:
#
# - [1] Forets, Marcelo, and Amaury Pouly. "Explicit error bounds for Carleman linearization."
#       arXiv preprint arXiv:1711.02552 (2017).
#
# - [2] Liu, J. P., Kolden, H. Ø., Krovi, H. K., Loureiro, N. F., Trivisa, K., & Childs,7
#       A. M. (2021). "Efficient quantum algorithm for dissipative nonlinear differential equations."
#       Proceedings of the National Academy of Sciences, 118(35). arXiv preprint arXiv:2011.03185.
#
# =========================================================================================================================

# --- Error bounds using a priori estimate from [1] ---
# These bounds use the supremum norm (p = Inf).

# see Theorem 4.2 in [1]; this bound is based on an a priori estimate
# of the norm of the exact solution x(t)
function error_bound_apriori(α, F₁, F₂; N)
    nF₂ = opnorm(F₂, Inf)
    μF₁ = logarithmic_norm(F₁, Inf)

    β = α * nF₂ / μF₁
    ε = t -> α * β^N * (exp(μF₁ * t) - 1)^N
    return ε
end

# see Theorem 4.2 in [1]
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

# --- Error bounds using power series method from [1] ---

# see Theorem 4.3 in [1], which uses the power series method
function error_bound_pseries(x₀, F₁, F₂; N)
    nx₀ = norm(x₀, Inf)
    nF₁ = opnorm(F₁, Inf)
    nF₂ = opnorm(F₂, Inf)
    β₀ = nx₀ * nF₂ / nF₁

    ε = t -> nx₀ * exp(nF₁ * t) / (1 - β₀ * (exp(nF₁ * t) - 1)) * (β₀ * (exp(nF₁ * t) - 1))^N
    return ε
end

# see Theorem 4.3 in [1]
function convergence_radius_pseries(x₀, F₁, F₂)
    nx₀ = norm(x₀, Inf)
    nF₁ = opnorm(F₁, Inf)
    nF₂ = opnorm(F₂, Inf)
    β₀ = nx₀ * nF₂ / nF₁

    T = (1 / nF₁) * log(1 + 1 / β₀)
    return T
end

# --- Error bounds using spectral abscissa from [2] ---
# These bounds use the spectral norm (p = 2).

# compute eigenvalues and sort them by increasing real part
function _error_bound_specabs_Re_λ₁(F₁; check=true)
    λ = eigvals(F₁; sortby=real)
    λ₁ = last(λ)
    Re_λ₁ = real(λ₁)
    if check && Re_λ₁ > 0
        throw(ArgumentError("expected Re(λ₁) ≤ 0, got $Re_λ₁"))
    end
    return Re_λ₁
end

# see Equation (2.2) in [2]
function _error_bound_specabs_R(x₀, F₂, Re_λ₁)
    nx₀ = norm(x₀, 2)
    nF₂ = opnorm(F₂, 2)
    R = nx₀ * nF₂ / abs(Re_λ₁)
    return R
end

# see Corollary 1 in [2]
function error_bound_specabs(x₀, F₁, F₂; N, check=true)
    Re_λ₁ = _error_bound_specabs_Re_λ₁(F₁; check=check)
    R = _error_bound_specabs_R(x₀, F₂, Re_λ₁)
    if check && R >= 1
        throw(ArgumentError("expected R < 1, got R = $R; try scaling the ODE"))
    end

    nx₀ = norm(x₀, 2)
    if iszero(Re_λ₁)
        nF₂ = opnorm(F₂, 2)
        ε = t -> nx₀ * (nx₀ * nF₂ * t)^N
    else
        ε = t -> nx₀ * R^N * (1 - exp(Re_λ₁ * t))^N  # Equation (4.29)
    end
    return ε
end

# see Corollary 1 in [2]
function convergence_radius_specabs(x₀, F₁, F₂; check=true)
    Re_λ₁ = _error_bound_specabs_Re_λ₁(F₁; check=check)

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
