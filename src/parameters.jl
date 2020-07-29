"""
            Objects for SPH mapping.


    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2020-07-03

"""

# using Distributed



"""
    struct mappingParameters
        x_lim::Vector{Float64}
        y_lim::Vector{Float64}
        z_lim::Vector{Float64}
        center::Vector{Float64}
        pixelSideLength::Float64
        pixelArea::Float64
        Npixels::Vector{Int64}
        x_size::Float64
        y_size::Float64
        z_size::Float64
    end

# Constructor:

    mappingParameters(;x_lim::Vector{Float64}   = [-1.0, -1.0],
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
"""
struct mappingParameters

    x_lim::Vector{Float64}
    y_lim::Vector{Float64}
    z_lim::Vector{Float64}
    center::Vector{Float64}
    pixelSideLength::Float64
    pixelArea::Float64
    Npixels::Vector{Int64}
    x_size::Float64
    y_size::Float64
    z_size::Float64
    boxsize::Float64
    periodic::Bool

    function mappingParameters(;x_lim::Vector{Float64}   = [-1.0, -1.0],
                                y_lim::Vector{Float64}   = [-1.0, -1.0],
                                z_lim::Vector{Float64}   = [-1.0, -1.0],
                                center::Vector{Float64}  = [-1.0, -1.0, -1.0],
                                x_size::Float64          =  -1.0,
                                y_size::Float64          =  -1.0,
                                z_size::Float64          =  -1.0,
                                pixelSideLength::Float64 =  -1.0,
                                Npixels::Int64           =   0,
                                boxsize::Float64         =  -1.0)


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

        # find the maximum extent of the map
        max_size = max(x_size, y_size)

        if (pixelSideLength == -1.0) & (Npixels != 0)
            pixelSideLength = max_size/Npixels
        elseif (pixelSideLength != -1.0) & (Npixels == 0)
            Npixels = floor(Int64, max_size/pixelSideLength)
            # recalculate pixelSideLenght to account for rounding
            pixelSideLength = max_size/Npixels
        else
            error("Please specify pixelSideLenght or number of pixels!")
        end

        pixelArea = pixelSideLength^2

        Npix = [ floor(x_size/pixelSideLength), 
                 floor(y_size/pixelSideLength),
                 floor(z_size/pixelSideLength) ]

        periodic = false 

        if boxsize != -1.0
            periodic = true
        end


        new(x_lim, y_lim, z_lim,
            center,
            pixelSideLength,
            pixelArea,
            Npix,
            x_size, y_size, z_size,
            boxsize, periodic)

    end
end
