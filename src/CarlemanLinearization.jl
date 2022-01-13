module CarlemanLinearization

# dependencies
include("init.jl")

# method
include("kronecker.jl")
include("linearization.jl")
include("error_bounds.jl")
include("solve.jl")

# exported methods and types
include("exports.jl")

end # module
