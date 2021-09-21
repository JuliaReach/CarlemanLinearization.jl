ENV["GKSwstype"] = "100"  # set 'GR environment' to 'no output' (for Travis CI)
using Documenter, CarlemanLinearization
using DynamicPolynomials

DocMeta.setdocmeta!(CarlemanLinearization, :DocTestSetup,
                   :(using CarlemanLinearization); recursive=true)

# generate Literate documentation
# include("generate.jl")

makedocs(
    format = Documenter.HTML(prettyurls = haskey(ENV, "GITHUB_ACTIONS"),  # disable for local builds
                             collapselevel = 1),
    sitename = "CarlemanLinearization.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "API Reference" => Any["Linearization" => "lib/linearization.md",
                               "Error bounds" => "lib/errors.md"],
        "References" => "references.md",
        "About" => "about.md"
    ]
)

# Deploy built documentation from Travis.
deploydocs(
    repo = "github.com/JuliaReach/CarlemanLinearization.jl.git",
    push_preview = true,
)
