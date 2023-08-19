"""
    synchrotron_SB_to_luminosity(map, pixelSideLength::Real)

Converts a map of synchrotron surface brightness `[erg/s/Hz/cm^2]` to synchrotron luminosity `[W/Hz]`.
Uses `pixelSideLength` as the diameter of a pixel in `[kpc]`.
"""
synchrotron_SB_to_luminosity(map::Matrix{<:Real}, pixelSideLength::Real) = 
    surface_brightness_to_luminosity(map, pixelSideLength, unit_factor=1.0e-7)


"""
    synchrotron_SB_to_luminosity(map, par::mappingParameters)

Converts a map of synchrotron surface brightness `[erg/s/Hz/cm^2]` to synchrotron luminosity `[W/Hz]`.
Uses `par` as the `mappingParameters` of the original map.
"""
synchrotron_SB_to_luminosity(map::Matrix{<:Real}, par::mappingParameters) =
         surface_brightness_to_luminosity(map, par, unit_factor=1.0e-7)


"""
    total_synch_luminosity_from_SB(map::Matrix{<:Real}, pixelSideLength::Real)

Computes the total synchrotron luminosity in `[W/Hz]` from a map of synchrotron surface brightness in `[erg/s/Hz/cm^2]`.
"""
function total_synch_luminosity_from_SB(map::Matrix{<:Real}, pixelSideLength::Real)
    # convert SB to Luminosity and return sum
    sum(synchrotron_SB_to_luminosity(map, pixelSideLength))
end

"""
    total_synch_luminosity_from_SB(map::Matrix{<:Real}, par::mappingParameters)

Computes the total synchrotron luminosity in `[W/Hz]` from a map of synchrotron surface brightness in `[erg/s/Hz/cm^2]`.
"""
total_synch_luminosity_from_SB(map::Matrix{<:Real}, par::mappingParameters) =
    total_synch_luminosity_from_SB(map, par.pixelSideLength)

"""
    total_synch_luminosity_from_SB(filename::String)

Computes the total synchrotron luminosity in `[W/Hz]` from a map of synchrotron surface brightness in `[erg/s/Hz/cm^2]`.
"""
function total_synch_luminosity_from_SB(filename::String)

    # read image and parameters
    im_jnu, par, snap_num, units = read_fits_image(filename)

    # convert SB to Luminosity
    return total_synch_luminosity_from_SB(im_jnu, par.pixelSideLength)
end