# Mapping SPH Data to Cartesian Grids

```@meta
CurrentModule = SPHtoGrid
DocTestSetup = quote
    using SPHtoGrid
end
```

Mapping SPH particles to a grid instead of simply plotting color-coded particle positions shows the actual gas quantities as they are used in an SPH code: weighted with a kernel, according to their distance to each other.
You can see this in the following plot, left are color-coded particle positions, right is the mean density in the SPH particles along the line of sight, interpolated to a grid.

![galaxy](assets/galaxy.png)

You can map SPH data to a grid using the function [`sphMapping`](@ref), which comes in two flavors: [CIC](@ref) and [TSC](@ref).

## Define parameters for mapping

To map the data you need to define the mapping parameters via the [`mappingParameters`](@ref) struct.

```@docs
mappingParameters
```

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

## CIC

For "Counts in Cell" (CIC) interpolation use the function `sphMapping` with these input values:

```@docs
sphMapping(::Array{<:Real}, ::Array{<:Real}, ::Array{<:Real}, ::Array{<:Real}, ::Array{<:Real})
```



### Select Kernel

You also need to choose the kernel you used in the simulation. For this you need to install the package [SPHKernels.jl](https://github.com/LudwigBoess/SPHKernels.jl). You can currently use these kernels:

```julia
k = Cubic()
k = Quintic()
k = WendlandC4()
k = WendlandC6()
k = WendlandC8()
```

Please see the SPHKernels [docs](https://ludwigboess.github.io/SPHKernels.jl/stable/) for more details.

### Mapping

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

Particles are mapped to a grid while also conserving the particle volume, following the algorithm described in [Dolag et. al. 2006](https://ui.adsabs.harvard.edu/link_gateway/2005MNRAS.363...29D/doi:10.1111/j.1365-2966.2005.09452.x).

### Weight functions

With the mapping you may decide to use a specific weighting function. For this you can pass the optional variable `Weights` in [`sphMapping`](@ref).

You can either use your own weight functions or use one of the built-in ones:

```@docs
part_weight_one
part_weight_physical
part_weight_emission
part_weight_spectroscopic
part_weight_XrayBand
```

### Units

As you have to handle unit conversion yourself please note that internally the image contains two components: 
`image` is the mapped pixel value and `wimage` is just the geometry weighting of the particle which will be applied if `reduce_image=true` in [`sphMapping`](@ref).

This means that the resulting images have the units:

` image = [Bin_Quant] * [Weights] * [pix^2] `

` wimage = [Weights] * [pix^2] `

### Helper Function

If you're lazy like me and don't want to go through the entire process of image reduction and saving the fits file by hand you can use

```@docs
map_it
```


## TSC

!!! note 
    Will be reworked for version 0.5!

