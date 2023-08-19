"""
    surface_brightness_to_luminosity(map::Matrix{<:Real}, pixelSideLength::Real; unit_factor::Real=1.0)

Converts a map of surface brightness to luminosity per pixel.
Uses `pixelSideLength` as the diameter of a pixel in `[kpc]`.
If `unit_factor` is provided it is multiplied to every pixel to perform unit conversion.
"""
function surface_brightness_to_luminosity(map::Matrix{<:Real}, pixelSideLength::Real; unit_factor::Real=1.0)

    # multiply by pixel area -> [erg/s]
    map_L = map .* (pixelSideLength * 3.085678e21)^2

    # convert to desired units
    map_L .*= unit_factor

    return map_L
end

"""
    surface_brightness_to_luminosity(map::Matrix{<:Real}, par::mappingParameters; unit_factor::Real=1.0)

Converts a map of surface brightness to luminosity per pixel.
If `unit_factor` is provided it is multiplied to every pixel to perform unit conversion.
"""
surface_brightness_to_luminosity(map::Matrix{<:Real}, par::mappingParameters; unit_factor::Real=1.0) =
    surface_brightness_to_luminosity(map, par.pixelSideLength; unit_factor)