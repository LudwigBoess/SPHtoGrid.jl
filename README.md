[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://LudwigBoess.github.io/SPHtoGrid.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://LudwigBoess.github.io/SPHtoGrid.jl/dev)
[![The MIT License](https://img.shields.io/badge/license-MIT-orange.svg)](LICENSE.md)
[![Build Status](https://travis-ci.org/LudwigBoess/SPHtoGrid.jl.svg?branch=master)](https://travis-ci.org/LudwigBoess/SPHtoGrid.jl)
[![codecov.io](https://codecov.io/gh/LudwigBoess/SPHtoGrid.jl/coverage.svg?branch=master)](https://codecov.io/gh/LudwigBoess/SPHtoGrid.jl?branch=master)

# SPHtoGrid.jl

This package maps SPH quantities to a cartesian grid.

You can map SPH data to a grid using the function `sphMapping`:

```julia
sphMapping( Pos, HSML, M, ρ, Bin_Quant;
            param::mappingParameters,
            kernel::SPHKernel,
            show_progress::Bool=true,
            conserve_quantities::Bool=false,
            parallel::Bool=false,
            filter_particles::Bool=true,
            dimensions::Int=2 )
```

### Setup

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

You also need to choose the kernel you used in the simulation. I implemented the following ones:

```julia
k = Cubic()
k = Quintic()
k = WendlandC4()
k = WendlandC6()
```

### Mapping

With the setup done you can now map (e.g.) density of your data using the function above as:

```julia
image = sphMapping(x, hsml, m, rho, rho, param=par, kernel=k)
```

Replacing the second `rho` with any other quantity would map that quantity of course.
Please note: This function doesn't do any unit conversion for you, so you need to convert to the desired units beforehand. See the chapter on unit conversion for usage.

Image now contains a 2D array with the binned data and can easily be plotted with `imshow()` from any plotting package of your choosing.

Per default the keyword `parallel = true` causes the run to use multiple processors. For this you need to start julia with `julia -p <N>` where `<N>` is the number of processors in your machine.

### Conserved quantities

With the latest release you can map the particles to a grid while also conserving the particle volume, following the algorithm described in [Dolag et. al. 2006](https://ui.adsabs.harvard.edu/link_gateway/2005MNRAS.363...29D/doi:10.1111/j.1365-2966.2005.09452.x).

This is switched off by default since it's slightly more expensive than simple mapping. If you want to use it simply call the mapping function with `conserve_quantities=true`.