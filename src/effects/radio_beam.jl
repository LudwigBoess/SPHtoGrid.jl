using Cosmology
import GadgetIO.AbstractGadgetHeader
using GadgetUnits
using Unitful

"""
    beam_in_kpc(θ_beam::Vector{Union{Real, Unitful.AbstractQuantity}}, 
                c::Cosmology.AbstractCosmology, z::Real)

Converts the beam from arcmin to kpc.
"""
function beam_in_kpc(θ_beam::Vector{<:Union{<:Real, <:Unitful.AbstractQuantity}}, 
                     c::Cosmology.AbstractCosmology, z::Real)

    smooth_size = zeros(2)
    for i = 1:2
        smooth_size[i] = 0.5 * ustrip(arcmin_to_kpc(c, θ_beam[i], z))
    end

    smooth_size
end


"""
    convert_Pnu_map_to_mJy_beam(map::Matrix{<:Real}, 
                                d_pixel::Real,
                                beam::Vector{Union{Real, Unitful.AbstractQuantity}}, 
                                c::Cosmology.AbstractCosmology, 
                                z::Real)

Converts a map from units [W / Hz] to [mJy / beam].

## Parameters:
- `map`: original map in [W / Hz].
- `d_pixel`: size of a pixel in ``kpc``.
- `beam`: dimensions of the beam in ``arcmin``.
- `c`: Cosmology used for conversion.
- `z`: Redshift of the image.
"""
function convert_Pnu_map_to_mJy_beam(map::Matrix{<:Real}, 
                                     d_pixel::Real,
                                     beam::Vector{<:Union{<:Real, <:Unitful.AbstractQuantity}}, 
                                     c::Cosmology.AbstractCosmology, 
                                     z::Real)

    if iszero(z)
        error("Conversion does not work for redshift zero!")
    end

    # convert the beam from arc-minutes to kpc
    kpc_beam = beam_in_kpc(beam, c, z)

    # area of beam ellipse in kpc
    beam_area = (π * (kpc_beam[1]*kpc_beam[2])) 

    # ratio between pixel and beam size
    pixel_factor = (d_pixel^2 / beam_area)

    # conversion step
    mJy_map = map ./ mJy_to_W(c, 1.0, z) ./ pixel_factor

    return mJy_map
end

"""
    convert_Pnu_map_to_mJy_beam(map::Matrix{<:Real}, 
                                d_pixel::Real,
                                beam::Union{Real, Unitful.AbstractQuantity}, 
                                c::Cosmology.AbstractCosmology, 
                                z::Real)

Converts a map from units [W / Hz] to [mJy / beam] for a circular beam.

## Parameters:
- `map`: original map in [W / Hz].
- `d_pixel`: size of a pixel in ``kpc``.
- `beam`: radius of the beam in ``arcmin``.
- `c`: Cosmology used for conversion.
- `z`: Redshift of the image.
"""
function convert_Pnu_map_to_mJy_beam(map::Matrix{<:Real}, 
                                     d_pixel::Real,
                                     beam::Union{<:Real, <:Unitful.AbstractQuantity},  
                                     c::Cosmology.AbstractCosmology, 
                                     z::Real)

    beam_arr = [beam, beam]
    convert_Pnu_map_to_mJy_beam(map, d_pixel, beam_arr, c, z)
end

"""
    convert_Pnu_map_to_mJy_beam(map::Matrix{<:Real}, 
                                d_pixel::Real,
                                beam::Union{T, Vector{T}}, 
                                h::AbstractGadgetHeader) where T::Union{Real, Unitful.AbstractQuantity}

Converts a map from units [W / Hz] to [mJy / beam] by using a Gadget header.

## Parameters:
- `map`: original map in [W / Hz].
- `d_pixel`: size of a pixel in ``kpc``.
- `beam`: radius/dimensions of the beam in ``arcmin``.
- `h`: Gadget header of simulation
"""
function convert_Pnu_map_to_mJy_beam(map::Matrix{<:Real}, 
                                     d_pixel::Real,
                                     beam::Union{T, Vector{T}}, 
                                     h::AbstractGadgetHeader) where T<:Union{<:Real, <:Unitful.AbstractQuantity}

    convert_Pnu_map_to_mJy_beam(map, d_pixel, beam, cosmology(h), h.z)
end