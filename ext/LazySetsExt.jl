module LazySetsExt

using CarlemanLinearization: generate_monomials
using LazySets: Hyperrectangle, dim, low, high
using IntervalArithmetic: interval

import CarlemanLinearization: lift_vector

"""
    lift_vector(X0, N)

Return a vector of monomials in `X0` (hyperrectangle) of degree at most `N`.
"""
function lift_vector(X0, N)
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
