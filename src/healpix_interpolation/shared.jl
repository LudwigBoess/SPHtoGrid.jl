function get_norm(_pos::Vector{Float64})

    Δx = 0.0
    @inbounds for dim  = 1:3
        Δx += _pos[dim]^2
    end

    return √(Δx)
end

function reduce_image_healpix(image_arr, Nside)

    image = reshape(image_arr[:,1], (Nside, 2Nside))
    w_image = reshape(image_arr[:,2], (Nside, 2Nside))

    for idx ∈ eachindex(image)
        if !isnan(w_image[idx]) && !iszero(w_image[idx])
            image[idx] /= w_image[idx]
        end
    end

    return image
end
