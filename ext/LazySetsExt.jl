module LazySetsExt

using CarlemanLinearization: generate_monomials
import CarlemanLinearization
@static if isdefined(Base, :get_extension)
    import LazySets
    using LazySets: Hyperrectangle, dim, low, high
    using IntervalArithmetic: inf, interval, sup
else
    import ..LazySets
    using ..LazySets: Hyperrectangle, dim, low, high
    using ..IntervalArithmetic: inf, interval, sup
end

"""
    lift_vector(X0, N)

Return a vector of monomials in `X0` (hyperrectangle) of degree at most `N`.
"""
function CarlemanLinearization.lift_vector(X0, N::Number)
    monoms = generate_monomials(dim(X0), N)
    nonfirst_monoms = Iterators.peel(monoms)[2]
    result = []
    intervals = [interval(low(X0, i), high(X0, i)) for i in 1:dim(X0)]
    for m in nonfirst_monoms
        push!(result, prod(intervals .^ m))
    end
    return Hyperrectangle(; low=[inf(i) for i in result], high=[sup(i) for i in result])
end

end  # module
