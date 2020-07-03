"""
            Kernels, Kernelvalues and -derivatives for sph mapping.


    Author: Ludwig Böss
    Contact: lboess@usm.lmu.de
    Created: 2018-12-12

"""

#using Distributed

abstract type SPHKernel end

"""
            Cubic
"""
# struct
struct Cubic <: SPHKernel
    n_neighbours::Int64
    norm_2D::Float64
    norm_3D::Float64
    function Cubic(n_neighbours::Int64=64)
        new(n_neighbours, 8.0/π, 8.0/π)
    end
end

# kernel values
@inline function kernel_value_2D(kernel::Cubic, u::Float64, h_inv::Float64)

    n = kernel.norm_3D * h_inv^2

    if u < 0.5
        return ( 1.0 - 6.0 * (1.0 - u) * u^2) * n
    elseif u < 1.0
        return ( 2.0 * (1.0 - u) * (1.0 - u) * (1.0 - u)) * n
    else
        return 0.
    end

end

@inline function kernel_value_3D(kernel::Cubic, u::Float64, h_inv::Float64)

    n = kernel.norm_3D * h_inv^3

    if u < 0.5
        return ( 1.0 + 6.0 * (u - 1.0) * u^2) * n
    elseif u < 1.0
        return ( 2.0 * (1.0 - u) * (1.0 - u) * (1.0 - u)) * n
    else
        return 0.
    end

end


"""
            Quintic
"""
# struct
struct Quintic <: SPHKernel
    n_neighbours::Int64
    norm_2D::Float64
    norm_3D::Float64
    function Quintic(n_neighbours::Int64=200)
        new(n_neighbours, 15309.0/(478.0*π), 2187.0/(40.0*pi))
    end
end

# kernel values
@inline function kernel_value_2D(kernel::Quintic, u::Float64, h_inv::Float64)

    n = kernel.norm_2D * h_inv^2

    if u < 1.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        u_m13 = ( 1.0/3.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 6.0 * u_m23*u_m23*u_m23*u_m23*u_m23 + 15.0 * u_m13*u_m13*u_m13*u_m13*u_m13 ) * n
    elseif u < 2.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 6.0 * u_m23*u_m23*u_m23*u_m23*u_m23 ) * n
    elseif u < 1.0
        u_m1  = ( 1.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 ) * n
    else
        return 0.
    end

end

@inline function kernel_value_3D(kernel::Quintic, u::Float64, h_inv::Float64)

    n = kernel.norm_3D * h_inv^3

    if u < 1.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        u_m13 = ( 1.0/3.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 6.0 * u_m23*u_m23*u_m23*u_m23*u_m23 + 15.0 * u_m13*u_m13*u_m13*u_m13*u_m13 ) * n
    elseif u < 2.0/3.0
        u_m1  = ( 1.0 - u )
        u_m23 = ( 2.0/3.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 - 6.0 * u_m23*u_m23*u_m23*u_m23*u_m23 ) * n
    elseif u < 1.0
        u_m1  = ( 1.0 - u )
        return ( u_m1*u_m1*u_m1*u_m1*u_m1 ) * n
    else
        return 0.
    end

end


"""
            Wendland C4
"""
# struct
struct WendlandC4 <: SPHKernel
    n_neighbours::Int64
    norm_2D::Float64
    norm_3D::Float64
    function WendlandC4(n_neighbours::Int64=200)
        new(n_neighbours,9.0/π, 495.0/(32.0 * π))
    end
end

# kernel values
@inline function kernel_value_2D(kernel::WendlandC4, u::Float64, h_inv::Float64)

    if u < 1.0
        n = kernel.norm_2D * h_inv^2
        u_m1 = 1.0 - u
        u_m1_2 = u_m1 * u_m1  # (1.0 - u)^2
        u_m1_4 = u_m1 * u_m1  # (1.0 - u)^4
        return ( u_m1_2*u_m1_4 * ( 1.0 + 6u + 35.0/3.0 * u^2 ) ) * n
    else
        return 0.
    end

end

@inline function kernel_value_3D(kernel::WendlandC4, u::Float64, h_inv::Float64)

    if u < 1.0
        n = kernel.norm_3D * h_inv^3
        u_m1 = 1.0 - u
        u_m1_2 = u_m1 * u_m1  # (1.0 - u)^2
        u_m1_4 = u_m1 * u_m1  # (1.0 - u)^4
        return ( u_m1_2*u_m1_4 * ( 1.0 + 6u + 35.0/3.0 * u^2 ) ) * n
    else
        return 0.
    end

end


"""
            Wendland C6
"""
# struct
struct WendlandC6 <: SPHKernel
    n_neighbours::Int64
    norm_2D::Float64
    norm_3D::Float64
    function WendlandC6(n_neighbours::Int64=295)
        new(n_neighbours, 78.0/(7.0*π), 1365.0/(64.0*π))
    end
end

# kernel values
@inline function kernel_value_2D(kernel::WendlandC6, u::Float64, h_inv::Float64)

    if u < 1.0
        n = kernel.norm_2D * h_inv^2
        u_m1 = 1.0 - u
        u_m1 = u_m1 * u_m1  # (1.0 - u)^2
        u_m1 = u_m1 * u_m1  # (1.0 - u)^4
        u_m1 = u_m1 * u_m1  # (1.0 - u)^8
        u2 = u*u
        return ( u_m1 * ( 1.0 + 8u + 25u2 + 32u2*u )) * n
    else
       return 0.0
   end

end

@inline function kernel_value_3D(kernel::WendlandC6, u::Float64, h_inv::Float64)

    if u < 1.0
        n = kernel.norm_3D * h_inv^3
        u_m1 = 1.0 - u
        u_m1 = u_m1 * u_m1  # (1.0 - u)^2
        u_m1 = u_m1 * u_m1  # (1.0 - u)^4
        u_m1 = u_m1 * u_m1  # (1.0 - u)^8
        u2 = u*u
        return ( u_m1 * ( 1.0 + 8u + 25u2 + 32u2*u )) * n
    else
       return 0.0
   end

end


# @inline function kernel_deriv(kernel::Cubic, u::Float64, h::Float64)
#
#     norm = 8.0/π
#     n = norm/h^4
#
#     if u < 0.5
#         return ( u * (18.0 * u - 12.0 )) * n
#     elseif u < 1.0
#         return ( -6.0 * (1.0 - u) * (1.0 - u) ) * n
#     else
#         return 0.
#     end
#
# end




# @inline function kernel_deriv(kernel::Quintic, u::Float64, h::Float64)
#
#     norm = ( 2187.0 / ( 40. * π))
#     n = norm/h^4
#
#     if u < 1.0/3.0
#         return ( -5.0 * ( 1.0 - u )^4 + 30.0 * ( 2.0/3.0 - u )^4  - 75.0 * ( 1.0/3.0 - u )^4 ) * n
#     elseif u < 2.0/3.0
#         return ( -5.0 * ( 1.0 - u )^4 + 30.0 * ( 2.0/3.0 - u )^4 - 75.0 ) * n
#     elseif u < 1.0
#         return ( -5.0 * ( 1.0 - u )^4 ) * n
#     else
#         return 0.
#     end
#
# end



# @inline function kernel_deriv(kernel::WendlandC4, u, h)
#
#     norm = 495.0/(32. * π)
#     n = norm/h^4
#
#     if u < 1.0
#         return ( -288.0/3.0 * ( 1. - u )^5 * u^2 - 56.0/3.0 * u * ( 1. - u)^5 ) * n
#     else
#         return 0.
#     end
#
# end



# @inline function kernel_deriv(kernel::WendlandC6, u::Float64, h::Float64)
#
#     norm = 1365.0/(64.0*π)
#     n = norm/h^4
#
#     if u < 1.0
#         return ( -22. * (1.0 - u)^7 * u * ( 16. * u^2 + 7. * u + 1. )) * n
#     else
#         return 0.
#     end
#
# end
