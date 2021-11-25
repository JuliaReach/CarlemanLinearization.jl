using Documenter, CarlemanLinearization
using DynamicPolynomials

DocMeta.setdocmeta!(CarlemanLinearization, :DocTestSetup,
                   :(using CarlemanLinearization); recursive=true)

makedocs(
    format = Documenter.HTML(prettyurls = haskey(ENV, "GITHUB_ACTIONS"),
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

deploydocs(
    repo = "github.com/JuliaReach/CarlemanLinearization.jl.git"
)
