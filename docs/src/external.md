# External Programs

SPHtoGrid.jl provides helper function for two external sph mapping Codes: Smac and P-Smac2.

## Smac

Smac is a SPH mapping Code by Klaus Dolag and others. The implementation is described in [Dolag et al. 2005](https://ui.adsabs.harvard.edu/link_gateway/2005MNRAS.363...29D/doi:10.1111/j.1365-2966.2005.09452.x).

Smac isn't public unfortunately. So these functions are mainly for my personal use.
If you do have access to Smac, here's a reference to what you can do.

SPHtoGrid.jl provides some functions to read the FITS and binary output of Smac.
To get the binary format you need to set `FILE_FORMAT = 1` in the parameter file.

### Writing the parameter file

You can write a Smac1 parameter file with

```@docs
write_smac1_par
```

### Reading image information

If you set `FILE_HEADER = 1` in the Smac parameter file you can read the information of the image header into a `Smac1ImageInfo` with:
```@docs
read_smac1_binary_image
```

For a Smac1 FITS file you can read the header into a [`mappingParameters`](@ref) struct with

```@docs
read_smac1_fits_info
```

### Reading the image

The image itself can be read with

```@docs
read_smac1_fits_image
read_smac1_binary_image
```

## P-Smac2

[P-Smac2](https://github.com/jdonnert/Smac2) by Julius Donnert is an advanced mapping code for a multitude of different quantities.

### Writing the parameter file

To run a mapping and plotting loop from a Julia script you need to update the parameter files on the fly.
The function `write_smac2_par` provides this functionality.

```@docs
write_smac2_par
```


### Reading image information
For a Smac1 FITS file you can read the header into a [`mappingParameters`](@ref) struct with

```@docs
read_smac2_info
```

### Reading the image

The image itself can be read with

```@docs
read_smac2_image
```