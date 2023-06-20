"""
    polarisation_fraction(Q_image, U_image, Iν_image, Iν_cutoff= 0.0)

Compute the polarisation fraction image from Stokes Q and U images.
"""
function polarisation_fraction(Q_image, U_image, Iν_image, Iν_cutoff= 0.0)
    
    # compute polarisation fraction
    Π_image = @. √(Q_image^2 + U_image^2) / Iν_image
    
    # set all pixels below cutoff to 0
    Π_image[Iν_image .<= Iν_cutoff] .= 0.0
    
    # return resulting image
    return Π_image
end

"""
    polarisation_angle(Q_image, U_image, Iν_image=nothing, Iν_cutoff= 0.0)

Compute the polarisation fraction image from Stokes Q and U images.
"""
function polarisation_angle(Q_image, U_image, Iν_image=nothing, Iν_cutoff= 0.0)

    # compute polarisation angle in degree
    ψ_image = @. 0.5atan(U_image / Q_image) |> rad2deg

    # if defined, set intensity cutoff
    if !isnothing(Iν_image)
        ψ_image[Iν_image.<=Iν_cutoff] .= 0.0
    end

    # return resulting image
    return ψ_image
end