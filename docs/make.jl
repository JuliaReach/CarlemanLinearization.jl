using Documenter, CarlemanLinearization, DocumenterCitations
using DynamicPolynomials

DocMeta.setdocmeta!(CarlemanLinearization, :DocTestSetup,
                    :(using CarlemanLinearization); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:alpha)

makedocs(; sitename="CarlemanLinearization.jl",
         format=Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true",
                                assets=["assets/aligned.css", "assets/citations.css"]),
         pagesonly=true,
         plugins=[bib],
         pages=["Home" => "index.md",
                "Bibliography" => "bibliography.md",
                "About" => "about.md"])

deploydocs(; repo="github.com/JuliaReach/CarlemanLinearization.jl.git",
           push_preview=true)
