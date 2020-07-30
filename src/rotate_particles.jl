using Base.Threads

"""
    euler_matrix(α::Real, β::Real, γ::Real)

Returns the rotation matrix A based on rotation along the euler angles α, β and γ corresponding to rotations around the x, y, and z-axis respectively.
"""
function euler_matrix(α::Real, β::Real, γ::Real)
    [
        (cos(β)*cos(γ))                         (-cos(β)*sin(γ))                        sin(β)
        (cos(α)*sin(γ) + cos(γ)*sin(α)*sin(β))  (cos(α)*cos(γ) - sin(α)*sin(β)*sin(γ))  (-cos(β)*sin(α))
        (sin(α)*sin(γ) - cos(α)*cos(γ)*sin(β))  (cos(γ)*sin(β) + cos(α)*sin(β)*sin(γ))  (cos(α)*cos(β))

    ]
end


"""
    rotate_3D_quantity(x, α, β, γ)

Rotates a 3D vector along the x-axis with angle α, y-axis with angle β and z-axis with angle γ.
"""
function rotate_3D_quantity(x::Array{<:Real}, α::Real, β::Real, γ::Real)

    A = euler_matrix(α, β, γ)

    return A * x
end


"""
    rotate_3D(x::AbstractArray, alpha::Real, beta::Real, gamma::Real)

Rotates and array of 3D positions around the euler angles α, β and γ corresponding to rotations around the x, y, and z-axis respectively.
α, β and γ need to be given in degrees.
"""
function rotate_3D(x::Array{<:Real}, alpha::Real, beta::Real, gamma::Real)

    α = deg2rad(alpha)
    β = deg2rad(beta)
    γ = deg2rad(gamma)

    N = length(x[:,1])
    ret = zeros(eltype(x[1,1]), N,3)
    @threads for i = 1:N
        @inbounds ret[i,:] = rotate_3D_quantity(x[i,:], α, β, γ)
    end

    return ret
end


"""
    rotate_3D!(x::AbstractArray, alpha::Real, beta::Real, gamma::Real)

Rotates and array of 3D positions around the euler angles α, β and γ corresponding to rotations around the x, y, and z-axis respectively.
α, β and γ need to be given in degrees.
"""
function rotate_3D!(x::Array{<:Real}, alpha::Real, beta::Real, gamma::Real)

    α = deg2rad(alpha)
    β = deg2rad(beta)
    γ = deg2rad(gamma)

    @threads for i = 1:length(x[:,1])
        @inbounds x[i,:] = rotate_3D_quantity(x[i,:], α, β, γ)
    end

    return x
end


"""
    rotate_to_xz_plane(x::AbstractArray)

Rotates an array of 3D positions into the xz-plane.
"""
rotate_to_xz_plane(x::Array{<:Real}) = [ x[:,1] x[:,3] x[:,2] ]

"""
    rotate_to_yz_plane(x::AbstractArray)

Rotates an array of 3D positions into the yz-plane.
"""
rotate_to_yz_plane(x::Array{<:Real}) = [ x[:,2] x[:,3] x[:,1] ]


"""
    project_along_axis(x::AbstractArray, projection_axis::Integer=3)

Projects and array of 3D along one of the principle axes.
projection_axis ∈ {1, 2, 3} => x, y, z axis.
"""
function project_along_axis(x::Array{<:Real}, projection_axis::Integer=3)
   
    # rotation to xy-plane -> nothing is done
    if projection_axis == 3
        return x
    end

    # rotation to xz-plane
    if projection_axis == 2
        return rotate_to_xz_plane(x)
    end

    # rotation to yz-plane
    if projection_axis == 1
        return rotate_to_yz_plane(x)
    end
end