using Documenter
using BaroclinicAdjustment

api = Any[
    "API" => "functions.md"
]

makedocs(
    sitename = "BaroclinicAdjustment",
    format = Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "API" => api,
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(repo = "github.com/simone-silvestri/BaroclinicAdjustment.jl.git", push_preview = true)
