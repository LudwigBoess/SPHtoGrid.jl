using Documenter
using SPHtoGrid

makedocs(
    sitename = "SPHtoGrid",
    format = Documenter.HTML(),
    modules = [SPHtoGrid]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/LudwigBoess/SPHtoGrid.jl.git"
)
