using Documenter
using BaroclinicAdjustment

closures = Any[
    "Horizontal closures" => "functions.md"
]

makedocs(
    sitename = "BaroclinicAdjustment",
    format = Documenter.HTML(),
    pages=[
        "Home" => "index.md",
        "Horizontal closures" => closures,
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(repo = "github.com/simone-silvestri/BaroclinicAdjustment.jl.git", push_preview = true)
