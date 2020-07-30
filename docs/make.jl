using Documenter
using SPHtoGrid

makedocs(
    sitename = "SPHtoGrid.jl",
    format = Documenter.HTML(),
    modules = [SPHtoGrid],
    pages = [
            "Table of Contents" => "index.md",
            "Install" => "install.md",
            "Mapping SPH Data" => "mapping.md",
            "Rotating Images" => "rotating.md",
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
