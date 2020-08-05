# Mapping SPH Data

```@meta
CurrentModule = SPHtoGrid
DocTestSetup = quote
    using SPHtoGrid
end
```

You can map SPH data to a grid using the function [`sphMapping`](@ref), which comes in two flavors: [CIC](@ref) and [TSC](@ref).


# Define parameters for mapping

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

If you are mapping a periodic box you also can give the keyword `boxsize` to enable periodic mapping.

```julia
par = mappingParameters(center=[x0, y0, z0], 
                        x_size=x_size, 
                        y_size=y_size,
                        z_size=z_size,
                        boxsize=boxsize,
                        Npixels=200)
```

# CIC

For "Cloud in a Cell" (CIC) interpolation use the function `sphMapping` with these input values:

```julia
function sphMapping(Pos::Array{<:Real}, HSML::Array{<:Real}, 
                    M::Array{<:Real}, ρ::Array{<:Real}, 
                    Bin_Quant::Array{<:Real},
                    Weights::Array{<:Real}=ρ;
                    param::mappingParameters,
                    kernel::SPHKernel [,
                    show_progress::Bool=true,
                    parallel::Bool=false,
                    filter_particles::Bool=true,
                    dimensions::Int=2])


    [...]

end
```

where `Pos` is the 3D positional data, `HSML` is the kernel support of each particle, `M` are the particle masses, `ρ` is the density of the particle, `Bin_Quant` is the qantity you want to bin and `Weights` is an array of weights. These weights can be calculated using one of the supplied [Weight functions](@ref). They default to density weighted.

## Select Kernel

You also need to choose the kernel you used in the simulation. For this you need to install the package [SPHKernels.jl](https://github.com/LudwigBoess/SPHKernels.jl). You can currently use these kernels:

```julia
k = Cubic()
k = Quintic()
k = WendlandC4()
k = WendlandC6()
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

With the latest release (v.0.3.0) the particles are mapped to a grid while also conserving the particle volume, following the algorithm described in [Dolag et. al. 2006](https://ui.adsabs.harvard.edu/link_gateway/2005MNRAS.363...29D/doi:10.1111/j.1365-2966.2005.09452.x).

## Weight functions

With the mapping you may decide to use a specific weighting function. For this you can pass the optional variable `Weights` in [`sphMapping`](@ref).

You can either use your own weight functions or use one of the built-in ones:

[`part_weight_one`](@ref) just returns an `Array` of ones.

[`part_weight_physical`](@ref) converts from pixel- to physical units.

[`part_weight_emission`](@ref) weights the contribution due to density and temperature of the particle.

[`part_weight_spectroscopic`](@ref) gives spectroscopic weighting, see Mazotta+ 04.

[`part_weight_XrayBand`](@ref) weights the particle due to its Xray emission in the defined energy band.


# TSC

A very simplistic approach to mapping SPH data is by using "Triangular Shaped Cloud" (TSC) interpolation. This has the advantage of only needing positional data to bin a quantity, so it works for mapping e.g. initial conditions that don't have values for `HSML` yet. To use this method call the [`sphMapping`](@ref) function with a reduced number of arguments:

```julia
sphMapping( Pos::Array{<:Real}, Bin_Quant::Array{<:Real};
            param::mappingParameters,
            show_progress::Bool=true,
            dimensions::Int=2)
```

Mapping then works the same as in the [CIC](@ref) case.

This only runs on a single core and will likely stay that way, as it is a "quick-and-dirty" method anyway.