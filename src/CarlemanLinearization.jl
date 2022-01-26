module CarlemanLinearization

# dependencies
include("init.jl")

# methods
include("kronecker.jl")
include("linearization.jl")
include("error_bounds.jl")
include("conversion.jl")

# exported methods and types
include("exports.jl")

end # module
