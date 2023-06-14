using Documenter, CarlemanLinearization
using DynamicPolynomials

DocMeta.setdocmeta!(CarlemanLinearization, :DocTestSetup,
                    :(using CarlemanLinearization); recursive=true)

makedocs(; sitename="CarlemanLinearization.jl",
         format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true",
                                assets=["assets/aligned.css"]),
         strict=true,
         pages=["Home" => "index.md",
                "API Reference" => Any["Linearization" => "lib/linearization.md",
                                       "Error bounds" => "lib/errors.md"],
                "References" => "references.md",
                "About" => "about.md"])

deploydocs(; repo="github.com/JuliaReach/CarlemanLinearization.jl.git",
           push_preview=true)
