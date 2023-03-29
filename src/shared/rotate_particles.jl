"""
    rotate_3D(x::Array{<:Real}, alpha::Real, beta::Real, gamma::Real)

Rotates and array of 3D positions around the euler angles α, β and γ corresponding to rotations around the x, y, and z-axis respectively.
α, β and γ need to be given in degrees.
"""
function rotate_3D(x::Array{<:Real}, alpha::Real, beta::Real, gamma::Real)

    # define rotation
    rot = RotXYZ(deg2rad(alpha), deg2rad(beta), deg2rad(gamma))

    return rot * x
end


"""
    rotate_3D!(x::Array{<:Real}, alpha::Real, beta::Real, gamma::Real)

Rotates and array of 3D positions around the euler angles α, β and γ corresponding to rotations around the x, y, and z-axis respectively.
α, β and γ need to be given in degrees.
"""
function rotate_3D!(x::Array{<:Real}, alpha::Real, beta::Real, gamma::Real)

    # define rotation
    rot = RotXYZ(deg2rad(alpha), deg2rad(beta), deg2rad(gamma))

    # rotate particles
    x = rot * x

    return x
end


"""
    rotate_to_xz_plane!(x::Array{<:Real})

Rotates an array of 3D positions into the xz-plane.
"""
function rotate_to_xz_plane!(x::Array{<:Real}) 

    @inbounds for i = 1:size(x,2)
        pos3   = copy(x[2,i])
        x[2,i] = x[3,i]
        x[3,i] = pos3
    end
    x
end


"""
    rotate_to_xz_plane!(x::Array{<:Real})

Rotates an array of 3D positions into the xz-plane.
"""
function rotate_to_xz_plane!(x::Array{<:Real}, x_in::Array{<:Real}) 

    @inbounds for i = 1:size(x,2)
        x[1,i] = x_in[1,i]
        x[2,i] = x_in[3,i]
        x[3,i] = x_in[2,i]
    end
    x
end

"""
    rotate_to_yz_plane(x::Array{<:Real})

Rotates an array of 3D positions into the yz-plane.
"""
function rotate_to_yz_plane!(x::Array{<:Real}) 

    @inbounds for i = 1:size(x,2)
        pos3   = copy(x[1,i])
        x[1,i] = x[2,i]
        x[2,i] = x[3,i]
        x[3,i] = pos3
    end
    x
end

"""
    rotate_to_yz_plane(x::Array{<:Real}, x_in::Array{<:Real})

Rotates an array of 3D positions into the yz-plane.
"""
function rotate_to_yz_plane!(x::Array{<:Real}, x_in::Array{<:Real}) 

    @inbounds for i = 1:size(x,1)
        x[1,i] = x_in[2,i]
        x[2,i] = x_in[3,i]
        x[3,i] = x_in[1,i]
    end
    x
end

"""
    project_along_axis!(x::Array{<:Real}, projection_axis::Integer=3)

Projects and array of 3D along one of the principle axes.
projection_axis ∈ {1, 2, 3} => x, y, z axis.
"""
function project_along_axis!(x::Array{<:Real}, projection_axis::Integer=3)
   
    # rotation to xy-plane -> nothing is done
    if projection_axis == 3
        return x
    end

    # rotation to xz-plane
    if projection_axis == 2
        return rotate_to_xz_plane!(x)
    end

    # rotation to yz-plane
    if projection_axis == 1
        return rotate_to_yz_plane!(x)
    end
end

"""
    project_along_axis!(x::Array{<:Real}, projection_axis::Integer=3)

Projects and array of 3D along one of the principle axes.
projection_axis ∈ {1, 2, 3} => x, y, z axis.
"""
function project_along_axis!(x::Array{<:Real}, x_in::Array{<:Real}, projection_axis::Integer=3)
   
    # rotation to xy-plane -> nothing is done
    if projection_axis == 3
        return x_in
    end

    # rotation to xz-plane
    if projection_axis == 2
        return rotate_to_xz_plane!(x, x_in)
    end

    # rotation to yz-plane
    if projection_axis == 1
        return rotate_to_yz_plane!(x, x_in)
    end
end

"""
    project_along_axis!(x::Array{<:Real}, projection_axis::Integer=3)

Projects and array of 3D along one of the principle axes.
projection_axis ∈ {1, 2, 3} => x, y, z axis.
"""
function project_along_axis(x::Array{<:Real}, projection_axis::Integer=3)
   
    # allocate new array
    ret = Array{eltype(x[1]),2}(undef, 3, size(x,2))

    # rotation to xy-plane -> nothing is done
    if projection_axis == 3
        return x
    end

    # rotation to xz-plane
    if projection_axis == 2
        return rotate_to_xz_plane!(ret, x)
    end

    # rotation to yz-plane
    if projection_axis == 1
        return rotate_to_yz_plane!(ret, x)
    end
end