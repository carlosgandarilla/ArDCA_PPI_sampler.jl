using Documenter
using ArDCA_PPI_sampler

makedocs(
    sitename = "ArDCA_PPI_sampler Documentation",
    modules  = [ArDCA_PPI_sampler],
    format   = Documenter.HTML(),
    pages    = [
        "Home" => "index.md",
        "Guide" => "guide.md",
        "API Reference" => "api.md",
    ]
)

deploydocs(
    repo = "github.com/carlosgandarilla/ArDCA_PPI_sampler.jl.git",
    branch = "gh-pages",
)
