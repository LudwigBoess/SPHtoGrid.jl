using GadgetIO, GadgetUnits
using Unitful, UnitfulAstro
using Cosmology

"""
    calculate_synchrotron_power( image_file::String, 
                                 θ₁::Real, θ₂::Real, 
                                 z::Real, c::Cosmology.AbstractCosmology,
                                 read_mode::Integer=1 )

Compute the synchrotron surface brightness in ``mJy`` and the total synchrotron power of an image 
Uses the geometric beam sizes `θ₁` and `θ₂` in arcseconds, the redshift `z` and a given cosmology `c` to do the cosmological calculations.

## Read modes:
* `1`: Read image of SPHtoGrid.jl mapping (default).
* `2`: Read Smac1 binary image.
* `3`: Read Smac2 FITS image.
"""
function calculate_synchrotron_power(image_file::String, 
                                     θ₁::Real, θ₂::Real, 
                                     z::Real, c::Cosmology.AbstractCosmology,
                                     read_mode::Integer=1)

    # read image from SPHtoGrid mapping
    if read_mode == 1
        map, par, snap_num, units = read_fits_image(image_file)
    
        # read Smac1 binary image
    elseif read_mode == 2
        map = read_smac1_binary_image(image_file)
        smac1_info = read_smac1_binary_info(image_file)
        smac1_center = [smac1_info.xcm, smac1_info.ycm, smac1_info.zcm] ./ 3.085678e21
        par = mappingParameters(center=smac1_center, 
                                x_size=smac1_info.boxsize_kpc,
                                y_size=smac1_info.boxsize_kpc,
                                z_size=smac1_info.boxsize_kpc,
                                Npixels=smac1_info.boxsize_pix)
    # read Smac2 FITS image
    elseif read_mode == 3
        map = read_smac2_image(image_file)
        par = read_smac2_info(image_file)
    # error handling
    else
        error("Read mode $read_mode not specified!")
    end

    # compute total intensity in erg / cm^2 / s / Hz
    # following Donnert&Brunetti 2014, Eq. A7
    I_ν =  (par.Npixels[1] * par.pixelSideLength)^2 / par.Npixels[1]^2 * sum(map)

    # convert to mJy 
    # 1 Jy = 10^{-23} erg / cm^2 / s / Hz
    I_ν *= 1.e26 

    # transform to observed frequency
    I_obs = I_ν / ( 1 + z)^3

    # convert beam size from arcseconds to arcminutes
    θ₁    = θ₁ * 1.0u"arcsecond" |> u"arcminute"
    θ₂    = θ₂ * 1.0u"arcsecond" |> u"arcminute"
    
    # physical beam diameter at the object
    θ₁_kpc = arcmin_to_kpc(c, θ₁, z) |> ustrip
    θ₂_kpc = arcmin_to_kpc(c, θ₂, z) |> ustrip

    # Surface brightness is Intensity * Beam area
    S_obs        = I_obs * ( π * 0.5θ₁_kpc * 0.5θ₂_kpc )
    
    # convert surface brightness in mJy to synchrotron power in W/Hz
    # P(ν) = S(ν) * 4π * dL^2
    # dL is the luminosity distance
    P_W          = mJy_to_W(c, S_obs, z)

    return S_obs, P_W
end