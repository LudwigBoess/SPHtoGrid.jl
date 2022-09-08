| **Documentation**                                                 | **Build Status**                                                                                | **Licence**                                                                                |
|:-----------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:| :-----------------------------------------------------------------------------------------------:|
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://LudwigBoess.github.io/SPHtoGrid.jl/stable) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LudwigBoess.github.io/SPHtoGrid.jl/dev) | [![Build Status](https://github.com/LudwigBoess/SPHtoGrid.jl/workflows/Run%20CI%20on%20master/badge.svg)](https://travis-ci.org/LudwigBoess/SPHtoGrid.jl) [![codecov.io](https://codecov.io/gh/LudwigBoess/SPHtoGrid.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/SPHtoGrid.jl?branch=master) |  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.md) |

# SPHtoGrid.jl

> :warning: **This package is in need of a significant overhaul and is not necessarily production ready! I use it mainly for my own work and some functionality is experimental! Be especially cautious with the effect functions! Also: The 3D mapping is currently broken and will (hopefully) be fixed soon!**: Please sanity-check everything before you use it in production!

This package maps SPH quantities to a cartesian grid. It is based on [Smac](https://ui.adsabs.harvard.edu/abs/2005MNRAS.363...29D/abstract) by [Klaus Dolag](https://www.usm.uni-muenchen.de/~dolag/) und [Smac2](https://github.com/jdonnert/Smac2) by Julius Donnert.

You can map SPH data to a grid using the function `sphMapping`:

```julia
function sphMapping( Pos::Array{<:Real}, HSML::Array{<:Real}, M::Array{<:Real}, 
                    Rho::Array{<:Real}, Bin_Quant::Array{<:Real}, 
                    Weights::Array{<:Real}=Rho;
                    param::mappingParameters,
                    kernel::AbstractSPHKernel,
                    show_progress::Bool=true,
                    parallel::Bool=false,
                    reduce_image::Bool=true,
                    filter_particles::Bool=true,
                    dimensions::Int=2,
                    calc_mean::Bool=false)


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

If you are mapping a periodic box you also can give the keyword `boxsize` to enable periodic mapping.

```julia
par = mappingParameters(center=[x0, y0, z0], 
                        x_size=x_size, 
                        y_size=y_size,
                        z_size=z_size,
                        boxsize=boxsize,
                        Npixels=200)
```

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

### Conserved quantities

With the latest release you can map the particles to a grid while also conserving the particle volume, following the algorithm described in [Dolag et. al. 2006](https://ui.adsabs.harvard.edu/link_gateway/2005MNRAS.363...29D/doi:10.1111/j.1365-2966.2005.09452.x).

## Weight functions

With the mapping you may decide to use a specivic weighting function. For this you can pass the optional variable `Weights` in `sphMapping`.

You can either use your own weight functions or use one of the built-in ones:

`part_weight_one` just returns an `Array` of ones.

`part_weight_physical` converts from pixel- to physical units.

`part_weight_emission` weights the contribution due to density and temperature of the particle.

`part_weight_spectroscopic` gives spectroscopic weighting, see Mazotta+ 04.

`part_weight_XrayBand` weights the particle due to its Xray emission in the defined energy band.