
"""
    get_quantities_splash( pos, weight, hsml, 
                       rho, m, len2pix::T) where T

Helper function to convert quantities to pixel units and the correct data type.
"""
function get_quantities_splash(pos, hsml, len2pix::T) where {T}

    hsml *= T(len2pix)
    hsml_inv = T(1.0) / hsml

    return T.(pos), hsml, hsml_inv
end