using FITSIO


"""
    write_fits_image(filename::String, image::Array{<:Real}, 
                            par::mappingParameters; 
                            units::String="[i.u.]",
                            snap::Integer=0)

Writes a mapped image to a FITS file and stored the essential mapping parameters in the header.
"""
function write_fits_image(filename::String, image::Array{<:Real},
    par::mappingParameters;
    units::String = "[i.u.]",
    snap::Integer = 0)

    z_slice_kpc = abs(par.z_lim[1]) + abs(par.z_lim[2])

    header_keys = ["SNAP",
        "CENTER_X", "CENTER_Y", "CENTER_Z",
        "XMIN", "XMAX", "YMIN", "YMAX", "ZMIN", "ZMAX",
        "BOXSIZE", "Z_SLICE",
        "NPIXELS", "PIX_SIZE",
        "UNITS"
    ]

    header_values = [snap,
        par.center[1], par.center[2], par.center[3],
        par.x_lim[1], par.x_lim[2],
        par.y_lim[1], par.y_lim[2],
        par.z_lim[1], par.z_lim[2],
        par.boxsize, z_slice_kpc,
        maximum(par.Npixels), par.pixelSideLength,
        units
    ]

    header_comments = ["snapshot number",
        "image center (x)", "image center (y)", "image center (z)",
        "x limit left", "x limit right",
        "y limit left", "y limit right",
        "z limit left", "z limit right",
        "size of the image", "depth of the integrated slice",
        "image resolution", "pixel size",
        "units of the image"]


    # construct header object
    header = FITSHeader(header_keys, header_values, header_comments)

    # write the FITS file
    f = FITS(filename, "w")

    write(f, image, header = header)


end

"""
    read_fits_image(filename::String)

Read a FITS file and return the image, mappingParameters and the snapshot number.

# Returns
- image: A 2D array with the image pixels 
- par:   mappingParameters used for the image 
- snap:  Number of the mapped snapshot
- units: A unit string of the image

# Example
image, par, snap_nr, unit_string = read_fits_image(filename)
"""
function read_fits_image(filename::String; verbose::Bool = false)

    if verbose
        @info "Reading image: $filename"
    end

    f = FITS(filename)

    if verbose
        @info "Opened file"
    end

    # read image
    image = read(f[1])

    if verbose
        @info "Read Image"
    end

    # read the header
    header = read_header(f[1])

    if verbose
        @info "Read Header"
    end

    # construct mappingParameters by now
    par = mappingParameters(x_lim = [header["XMIN"], header["XMAX"]],
        y_lim = [header["YMIN"], header["YMAX"]],
        z_lim = [header["ZMIN"], header["ZMAX"]],
        Npixels = Int64(header["NPIXELS"]),
        boxsize = header["BOXSIZE"]
    )

    return image, par, header["SNAP"], header["UNITS"]
end