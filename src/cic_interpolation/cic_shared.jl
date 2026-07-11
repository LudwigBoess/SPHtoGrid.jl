"""
    get_weight_per_pixel(distr_area, distr_weight,
                                n_tot_pix, n_distr_pix,
                                dxdy, u, hsml_inv,
                                kernel)

Computes the area/volume contribution for each pixel and evaluates the kernel at the pixel center.
"""
function get_weight_per_pixel(distr_area, distr_weight,
                                n_tot_pix, n_distr_pix,
                                dxdy, u, hsml_inv,
                                kernel)

    # store are of pixel the particle contributes to
    _A = dxdy
    distr_area += dxdy

    # count up total pixels 
    n_tot_pix += 1

    # if pixel center is within kernel
    if u <= 1
        # evaluate kernel
        _wk = 𝒲(kernel, u, hsml_inv)

        # count up distributed weigh
        distr_weight += _wk * dxdy

        # if pixel center is within the kernel
        # count that pixel contributes to pixel
        n_distr_pix += 1
    else
        _wk = 0.0
    end

    return _A, _wk,
    distr_area, distr_weight,
    n_tot_pix, n_distr_pix
end

"""
    function pix_index_min_max(x::T, hsml::T, n_pixels::Int64) where T

Calculates the minimum and maximum pixel to which a particle contributes.
"""
function pix_index_min_max(x::T, hsml::T, n_pixels::Int64) where T

    iMin = max(floor(Integer, x - hsml), 0)
    iMax = min(floor(Integer, x + hsml), n_pixels-1)

    return iMin, iMax
end


"""
    function get_dxyz(x::Real, hsml::Real, i::Integer)

Calculates the extent of the current particle size in units of pixels.
"""
function get_dxyz(x::Real, hsml::Real, i::Integer)
    min(x + hsml, i + 1) - max(x - hsml, i)
end


"""
    function get_x_dx(x, hsml, i)


"""
@inline function get_x_dx(x, hsml, i)

    dx = get_dxyz(x, hsml, i)
    x_dist = x - i - 0.5

    return x_dist, dx
end


"""
    function get_xyz( pos, hsml, k::Integer,
                      par::mappingParameters)

Calculates `x, y, z` position in units of pixels and performs periodic mapping, if required.
"""
function get_xyz( pos::Array{T},
                  par::mappingParameters) where T

    x, y, z = pos
    # convert from code units to pixels
    x *= par.len2pix
    y *= par.len2pix
    z *= par.len2pix

    x += 0.5 * par.Npixels[1]
    y += 0.5 * par.Npixels[2]
    z += 0.5 * par.Npixels[3]
    
    return x, y, z

end



"""
    function update_image( image::Real, w_image::Real, 
                           wk::Real, bin_q::Real, 
                           geometry_norm::Real )

Applies the different contributions to the image and the weight image.
"""
function update_image!(image::Array{Float64}, idx::Integer,
                        pix_weight::Float64,
                        bin_q::Union{Float64, Vector{Float64}} )

    image[idx, end] += pix_weight
    for i = 1:length(bin_q)
        image[idx,i] += bin_q[i] * pix_weight
    end
    
    return image
end

"""
    faraday_rotate_pixel!(image::Array{Float64}, idx::Integer, 
                          pRM::Float64, pix_weight::Float64)

Faraday rotates the current pixel state based on the contribution from particle `p`.
"""
function faraday_rotate_pixel!(image::Array{Float64}, idx::Integer, 
                                pRM::Float64, pix_weight::Float64,
                                stokes::Bool)

    # calculate contribution of particle RM to current pixel
    _RM = pRM * pix_weight

    # reduce to minimum angle
    _RM = mod(_RM, π)

    # reconstruct current polarisation angle
    # if Stokes parameters are supposed to be mapped
    if stokes
        # convention: 1: stokes_Q 
        #             2: stokes_U
        Q = image[idx, 1]
        U = image[idx, 2]

        # construct total polarized emission in pixel
        Ipol = √(Q^2 + U^2)

        # construct polarisation angle in pixel 
        ψ = 0.5atan(U / Q)

        # save rotated pixel
        image[idx, 1] = Ipol * cos(2(ψ + _RM))
        image[idx, 2] = Ipol * sin(2(ψ + _RM))
    end
    
    return image
end


function free_memory(x, hsml, m, rho, bin_q, weights)
    x = hsml = m = rho = bin_q = weights = nothing
    GC.gc()
end

"""
    mass_conservation_report(image, M, Rho, Weights, len2pix; ndim=2)

Logs how much of the mappable quantity landed on the grid, as a
grid-vs-particle comparison. The CIC deposition sums the weight image to a
closed form per particle:

    Σ_pixels pix_weight = len2pix^(ndim+1) · Weights·M/Rho

(`ndim=2` → `len2pix³` for a 2D map, `ndim=3` → `len2pix⁴` for a 3D cube), so
the total of the weight image (`image[:, end]`) divided by that factor
recovers `Σ Weights·M/Rho` over the particles actually deposited. For the
default density weighting (`Weights == Rho`) that is the total mass.

Because each kept particle is renormalised over its *in-grid* pixels, this
quantity is conserved to machine precision for every particle overlapping ≥1
in-grid pixel centre (edge tails are piled onto the near edge, not lost). A
ratio < 1 therefore measures *coverage* loss — particles off-grid, filtered
out, or too small to hit a pixel centre — and NOT interpolation/resolution
error. Returns `(on_grid, expected)` in code-mass units.
"""
function mass_conservation_report(image::AbstractArray, M, Rho, Weights, len2pix; ndim::Integer=2)

    on_grid  = sum(@view image[:, end]) / len2pix^(ndim+1)   # Σ Weights·M/Rho deposited
    expected = sum(Weights .* M ./ Rho)                      # Σ Weights·M/Rho of inputs
    ratio    = iszero(expected) ? 1.0 : on_grid / expected

    @info "Mass conservation (grid vs. particles):"
    @info "  On grid:   $(on_grid)"
    @info "  Expected:  $(expected)"
    @info "  Ratio:     $(ratio)"
    @info "  Rel. loss: $(1.0 - ratio)  (coverage loss, not interpolation error)"

    return on_grid, expected
end