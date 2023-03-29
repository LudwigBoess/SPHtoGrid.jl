"""
    get_d_hsml( dx::Real, dy::Real, hsml_inv::Real )

Computes the distance in 2D to the pixel center in units of the kernel support.
"""
function get_d_hsml(dx::T, dy::T, hsml_inv::T) where T
    √( dx*dx + dy*dy ) * hsml_inv
end

"""
    get_d_hsml_3D( dx::Real, dy::Real, dz::Real, hsml_inv::Real )

Computes the distance in 3D to the pixel center in units of the kernel support.
"""
function get_d_hsml(dx::T, dy::T, dz::T, hsml_inv::T) where {T}
    √(dx * dx + dy * dy + dz * dz) * hsml_inv
end