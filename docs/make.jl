using Documenter
using SPHtoGrid

makedocs(
    sitename = "SPHtoGrid.jl",
    doctest =false,
    format=Documenter.HTML(),
    modules = [SPHtoGrid],
    pages = [
            "Table of Contents" => "index.md",
            "Install" => "install.md",
            "Mapping SPH Data" => "mapping.md",
            "HealPix maps" => "healpix.md",
            "Examples" => "effects.md",
            "Rotating Images" => "rotating.md",
            "Saving/Loading Images" => "io.md",
            "External Programs" => "external.md",
            "API reference" => "api.md"
            ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/LudwigBoess/SPHtoGrid.jl.git"
)
