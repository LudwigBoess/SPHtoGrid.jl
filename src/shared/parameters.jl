"""
            Objects for SPH mapping.


    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2020-07-03

"""

# using Distributed


struct mappingParameters{T}

    x_lim::Vector{T}
    y_lim::Vector{T}
    z_lim::Vector{T}
    center::Vector{T}
    halfsize::Vector{T}
    len2pix::T
    pixelSideLength::T
    Npixels::Vector{Int64}
    boxsize::T
    periodic::Bool

end

"""
    mappingParameters( T::DataType=Float64;
                       x_lim::Vector{Float64}   = [-1.0, -1.0],
                       y_lim::Vector{Float64}   = [-1.0, -1.0],
                       z_lim::Vector{Float64}   = [-1.0, -1.0],
                       center::Vector{Float64}  = [-1.0, -1.0, -1.0],
                       x_size::Float64          =  -1.0,
                       y_size::Float64          =  -1.0,
                       z_size::Float64          =  -1.0,
                       pixelSideLength::Float64 =  -1.0,
                       Npixels::Int64           =   0)

Parameter object for sph to grid mapping. Define either `*_lim`, or `center` and `*_size`.
Resolution is defined by `pixelSideLength` or `Npixels`.

Non-square maps are supported: pixels are always square in physical units (one
`pixelSideLength` for all axes), so different extents in `x`, `y` and `z` simply
result in different pixel counts per axis (`Npixels` becomes a per-axis vector).
When `Npixels` is passed it sets the resolution along the largest in-plane
(x/y) extent and the remaining axes are scaled to keep pixels square.
"""
function mappingParameters( T::DataType=Float64;
                            x_lim::Vector{<:Real}   = [-1.0, -1.0],
                            y_lim::Vector{<:Real}   = [-1.0, -1.0],
                            z_lim::Vector{<:Real}   = [-1.0, -1.0],
                            center::Vector{<:Real}  = [-1.0, -1.0, -1.0],
                            x_size::Real            =  -1.0,
                            y_size::Real            =  -1.0,
                            z_size::Real            =  -1.0,
                            pixelSideLength::Real   =  -1.0,
                            Npixels::Integer        =   0,
                            boxsize::Real           =  -1.0)


    # calculate limits if center position and sizes are given
    if ( x_lim == [-1.0, -1.0] && ( y_lim == [-1.0, -1.0] && z_lim == [-1.0, -1.0] ))

        if (center != [-1.0, -1.0, -1.0]) && ( x_size != -1.0 && (y_size != -1.0 && z_size != -1.0) )

            x_lim = [ center[1] - 0.5x_size, center[1] + 0.5x_size ]
            y_lim = [ center[2] - 0.5y_size, center[2] + 0.5y_size ]
            z_lim = [ center[3] - 0.5z_size, center[3] + 0.5z_size ]

        else
            error("Giving a center position requires extent in x, y and z direction.")
        end
    end

    # calculate side lengths from limits
    if x_size == -1.0
        x_size = x_lim[2] - x_lim[1]
    end
    if y_size == -1.0
        y_size = y_lim[2] - y_lim[1]
    end
    if z_size == -1.0
        z_size = z_lim[2] - z_lim[1]
    end

    if center == [-1.0, -1.0, -1.0]
        center[1] = x_lim[1] + 0.5*x_size
        center[2] = y_lim[1] + 0.5*y_size
        center[3] = z_lim[1] + 0.5*z_size
    end

    # A single pixel side length is used for every axis, so pixels stay square
    # in physical units even for non-square maps — each axis then just gets its
    # own pixel count. When `Npixels` is given it sets the resolution along the
    # largest in-plane (x/y) extent; the other axes are scaled to match.
    if (pixelSideLength != -1.0) & (Npixels != 0)
        error("Please specify pixelSideLength or number of pixels!")
    elseif (pixelSideLength == -1.0) & (Npixels == 0)
        error("Please specify pixelSideLength or number of pixels!")
    end

    max_size = max(x_size, y_size)

    if pixelSideLength == -1.0
        pixelSideLength = max_size / Npixels
    end

    # number of pixels per axis from the uniform pixel size
    Npix = [ max(1, round(Int64, x_size / pixelSideLength)),
             max(1, round(Int64, y_size / pixelSideLength)),
             max(1, round(Int64, z_size / pixelSideLength)) ]

    periodic = false 

    if boxsize != -1.0
        periodic = true
    end

    halfsize = 0.5 .* [ x_size, y_size, z_size ]

    len2pix = 1.0/pixelSideLength

    mappingParameters{T}(x_lim, y_lim, z_lim,
                         center,
                         halfsize,
                         len2pix,
                         pixelSideLength,
                         Npix,
                         boxsize, periodic)

end
