# SPHtoGrid.jl Documenation

```@meta
CurrentModule = SPHtoGrid
DocTestSetup = quote
    using SPHtoGrid
end
```

```@contents
Pages = ["index.md"]
Depth = 3
```

## SPHtoGrid.jl


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

### Setup

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



External Programs
-----------------
GadJet.jl provides helper function for two external sph mapping Codes: Smac and P-Smac2.

### P-Smac2
P-Smac2 by Julius Donnert (https://github.com/jdonnert/Smac2) is an advanced mapping code for a multitude of different quantities. To run a mapping and plotting loop from a Julia script you need to update the parameter files on the fly.
The function `write_smac2_par` provides this functionality.

```julia
write_smac2_par(x, y, z,
                euler_angle_0, euler_angle_1, euler_angle_2,
                xy_size, z_depth, xy_pix::Int64,
                input_file, output_file, path,
                effect_module::Int64=0, effect_flag::Int64=0)
```

### Smac

Smac is a SPH mapping Code by Klaus Dolag and others. The implementation is described in Dolag et al. 2005 (https://ui.adsabs.harvard.edu/link_gateway/2005MNRAS.363...29D/doi:10.1111/j.1365-2966.2005.09452.x)

Smac isn't public unfortunately. So these functions are mainly for my personal use.
If you do have access to Smac, here's a reference to what you can do.

GadJet.jl provides some functions to read the binary output of Smac, as I personally prefer that over the FITS output.
To get the binary format you need to set `FILE_FORMAT = 1` in the parameter file.

#### Reading image information

If you set `FILE_HEADER = 1` in the Smac parameter file you can read the information of the image header into a `Smac1ImageInfo` object like so:

```julia
info = read_smac1_binary_info(filename)
```

The `Smac1ImageInfo` object contains the following information

```julia
struct Smac1ImageInfo

    snap::Int32                 # number of input snapshot
    z::Float32                  # redshift of snapshot
    m_vir::Float32              # virial mass of halo
    r_vir::Float32              # virial radius of halo
    xcm::Float32                # x coordinate of image center
    ycm::Float32                # y coordinate of image center
    zcm::Float32                # z coordinate of image center
    z_slice_kpc::Float32        # depth of the image in kpc
    boxsize_kpc::Float32        # xy-size of the image in kpc
    boxsize_pix::Float32        # xy-size of the image in pixels
    pixsize_kpc::Float32        # size of one pixel in kpc
    xlim::Array{Float64,1}      # x limits of image
    ylim::Array{Float64,1}      # y limits of image
    zlim::Array{Float64,1}      # z limits of image
    units::String               # unitstring of image

end

```

#### Reading the image

The image itself can be read with

```julia
image = read_smac1_binary_image(filename)
```

This will return an `Array{Float32,2}` with the pixel values. You can pass this to any imshow function of your favorite plotting package.


# Index
```@index
```

# API Reference

```@autodocs
Modules = [SPHtoGrid]
Order   = [:type, :function]
```