| **Documentation**                                                 | **Build Status**                                                                                | **Licence**                                                                                |**Citation**                                                                                |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:| :-----------------------------------------------------------------------------------------------:| :-----------------------------------------------------------------------------------------------:|
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://LudwigBoess.github.io/SPHtoGrid.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LudwigBoess.github.io/SPHtoGrid.jl/dev) | [![Build Status](https://github.com/LudwigBoess/SPHtoGrid.jl/workflows/Run%20CI%20on%20master/badge.svg)](https://travis-ci.org/LudwigBoess/SPHtoGrid.jl) [![codecov.io](https://codecov.io/gh/LudwigBoess/SPHtoGrid.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/SPHtoGrid.jl?branch=master) |  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.md) | [![DOI](https://zenodo.org/badge/276872260.svg)](https://zenodo.org/badge/latestdoi/276872260) |

# SPHtoGrid.jl

This package maps SPH quantities to a cartesian grid or a healpix sphere. It is based on [Smac](https://ui.adsabs.harvard.edu/abs/2005MNRAS.363...29D/abstract) by [Klaus Dolag](https://www.usm.uni-muenchen.de/~dolag/) und [Smac2](https://github.com/jdonnert/Smac2) by Julius Donnert.

Please see the [Documentation](https://LudwigBoess.github.io/SPHtoGrid.jl/dev) for details.

# Quickstart

You can map SPH data to a grid using the function `sphMapping`:

```julia
function sphMapping(Pos::Array{<:Real}, HSML::Array{<:Real}, M::Array{<:Real}, 
                    Rho::Array{<:Real}, Bin_Quant::Array{<:Real}, 
                    Weights::Array{<:Real}=Rho;
                    param::mappingParameters,
                    kernel::AbstractSPHKernel,
                    show_progress::Bool=true,
                    parallel::Bool=false,
                    reduce_image::Bool=true,
                    return_both_maps::Bool=false,
                    dimensions::Int=2,
                    calc_mean::Bool=false,
                    sort_z::Bool=false)


    [...]

end
```

## Define parameters for mapping

To map the data you need to define the mapping parameters via the `mappingParameters` object.
One way to set this up is by defining the limits of the map as

```julia
par = mappingParameters(xlim=[xmin, xmax],
                        ylim=[ymin, ymax],
                        zlim=[zmin, zmax],
                        Npixels=200)
```

or give a center position and the size in each direction

```julia
par = mappingParameters(center=[x0, y0, z0], 
                        x_size=x_size, 
                        y_size=y_size,
                        z_size=z_size,
                        Npixels=200)
```

Instead of Npixels you can also give the keyword argument `pixelSideLength` if you prefer to define your image that way.


## Select Kernel

You also need to choose the kernel you used in the simulation. For this you need to install the package [SPHKernels.jl](https://github.com/LudwigBoess/SPHKernels.jl). You can currently use these kernels:

```julia
k = Cubic()
k = Quintic()
k = WendlandC4()
k = WendlandC6()
k = WendlandC8()
```

Please see the SPHKernels [docs](https://ludwigboess.github.io/SPHKernels.jl/stable/) for more details.

## Mapping

With the setup done you can now map (e.g.) density of your data using the function above as:

```julia
image = sphMapping(x, hsml, m, rho, rho, param=par, kernel=k)
```

Replacing the second `rho` with any other quantity would map that quantity of course.
Please note: This function doesn't do any unit conversion for you, so you need to convert to the desired units beforehand. You can do this e.g. with [GadgetUnits.jl](https://github.com/LudwigBoess/GadgetUnits.jl).

Image now contains a 2D array with the binned data and can easily be plotted with `imshow()` from any plotting package of your choosing.

The keyword `parallel = true` causes the run to use multiple processors. For this you need to start julia with `julia -p <N>` where `<N>` is the number of processors in your machine, or define

```julia
using Distributed
addprocs(8)

# now you can load SPHtoGrid
using SPHtoGrid
```

## Conserved quantities

Particles are mapped to a grid while also conserving the particle volume, following the algorithm described in [Dolag et. al. 2006](https://ui.adsabs.harvard.edu/link_gateway/2005MNRAS.363...29D/doi:10.1111/j.1365-2966.2005.09452.x).
