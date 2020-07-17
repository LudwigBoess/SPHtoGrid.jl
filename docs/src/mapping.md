# Mapping SPH Data

```@meta
CurrentModule = SPHtoGrid
DocTestSetup = quote
    using SPHtoGrid
end
```

You can map SPH data to a grid using the function [`sphMapping`](@ref):

```julia
sphMapping( Pos, HSML, M, ρ, Bin_Quant;
            param::mappingParameters,
            kernel::SPHKernel,
            show_progress::Bool=true,
            conserve_quantities::Bool=false,
            parallel::Bool=false,
            filter_particles::Bool=true,
            dimensions::Int=2)
```

## Define parameters for mapping

To map the data you need to define the mapping parameters via the [`mappingParameters`](@ref) object.
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

You also need to choose the kernel you used in the simulation. I implemented the following ones:

```julia
k = Cubic()
k = Quintic()
k = WendlandC4()
k = WendlandC6()
```

Feel free to add new ones and issue pull requests!

## Mapping

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