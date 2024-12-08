module LazySetsExt

using CarlemanLinearization: generate_monomials
@static if isdefined(Base, :get_extension)
    import LazySets  # needed for `require`
    using LazySets: Hyperrectangle, dim, low, high
    using IntervalArithmetic: interval
else
    import ..LazySets
    using ..LazySets: Hyperrectangle, dim, low, high
    using ..IntervalArithmetic: interval
end

"""
    lift_vector(X0, N)

Return a vector of monomials in `X0` (hyperrectangle) of degree at most `N`.
"""
function _lift_vector(X0, N)
    monoms = generate_monomials(dim(X0), N)
    nonfirst_monoms = Iterators.peel(monoms)[2]
    result = []
    intervals = [interval(low(X0, i), high(X0, i)) for i in 1:dim(X0)]
    for m in nonfirst_monoms
        push!(result, prod(intervals .^ m))
    end
    return Hyperrectangle(; low=[i.lo for i in result], high=[i.hi for i in result])
end

end  # module
