module CarlemanLinearization

# dependencies
include("init.jl")

# method
include("linearization.jl")
include("error_bounds.jl")

# exported methods and types
include("exports.jl")

end # module
