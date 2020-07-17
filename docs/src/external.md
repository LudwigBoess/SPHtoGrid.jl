# External Programs

SPHtoGrid.jl provides helper function for two external sph mapping Codes: Smac and P-Smac2.

## P-Smac2

[P-Smac2](https://github.com/jdonnert/Smac2) by Julius Donnert is an advanced mapping code for a multitude of different quantities. To run a mapping and plotting loop from a Julia script you need to update the parameter files on the fly.
The function `write_smac2_par` provides this functionality.

```julia
write_smac2_par(x, y, z,
                euler_angle_0, euler_angle_1, euler_angle_2,
                xy_size, z_depth, xy_pix::Int64,
                input_file, output_file, path,
                effect_module::Int64=0, effect_flag::Int64=0)
```

## Smac

Smac is a SPH mapping Code by Klaus Dolag and others. The implementation is described in [Dolag et al. 2005](https://ui.adsabs.harvard.edu/link_gateway/2005MNRAS.363...29D/doi:10.1111/j.1365-2966.2005.09452.x).

Smac isn't public unfortunately. So these functions are mainly for my personal use.
If you do have access to Smac, here's a reference to what you can do.

SPHtoGrid.jl provides some functions to read the binary output of Smac, as I personally prefer that over the FITS output.
To get the binary format you need to set `FILE_FORMAT = 1` in the parameter file.

### Reading image information

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

### Reading the image

The image itself can be read with

```julia
image = read_smac1_binary_image(filename)
```

This will return an `Array{Float32,2}` with the pixel values. You can pass this to any imshow function of your favorite plotting package.