"""
    particle_in_image(x::Real, y::Real, z::Real, hsml::Real,
                                halfXsize::Real, halfYsize::Real, halfZsize::Real)

Checks if a periodically mapped particle is still in the image frame.
"""
@inline function particle_in_image( x::Real, y::Real, z::Real, hsml::Real,
                                    halfsize::Array{T},
                                    halfpixel::Real) where {T}

    if (( x - hsml + halfpixel ) > halfsize[1] || ( x + hsml - halfpixel ) < -halfsize[1] ||
        ( y - hsml + halfpixel ) > halfsize[2] || ( y + hsml - halfpixel ) < -halfsize[2] ||
        ( z - hsml + halfpixel ) > halfsize[3] || ( z + hsml - halfpixel ) < -halfsize[3] )

        return false
    else
        return true
    end
end

"""
    particle_in_image(x::Real, y::Real, z::Real)

Checks if a periodically mapped particle is still in the image frame.
"""
@inline function particle_in_image(x::Real, y::Real, z::Real, halfsize::Array{T}) where T

    if (x > halfsize[1] || x < -halfsize[1] ||
        y > halfsize[2] || y < -halfsize[2] ||
        z > halfsize[3] || z < -halfsize[3])
        return false
    else
        return true
    end
end

"""
    particle_in_image( pos::Array{T}, hsml::T,
                       halfsize::Array{Float64}) where T

Checks if a periodically mapped particle is still in the image frame.
"""
@inline function particle_in_image( pos::Vector{T}, hsml::T,
                                    halfsize::Vector{Float64},
                                    halfpixel) where T

    if ( ( pos[1] - hsml + halfpixel ) > halfsize[1] || ( pos[1] + hsml - halfpixel ) < -halfsize[1]  ||
         ( pos[2] - hsml + halfpixel ) > halfsize[2] || ( pos[2] + hsml - halfpixel ) < -halfsize[2]  ||
         ( pos[3] - hsml + halfpixel ) > halfsize[3] || ( pos[3] + hsml - halfpixel ) < -halfsize[3]  )

        return false
    else
        return true
    end
end

"""
    add_subtr_box(pos::T, boxsize::T) where T

Helper function to add or subtract the boxsize from the given position
"""
@inline function add_subtr_box(pos::T, boxsize::T) where T
    if pos > 0.0
        return pos - boxsize
    else
        return pos + boxsize
    end
end

"""
    find_position_periodic( pos, k, boxsize )

Performs a periodic mapping of the particle position.
"""
@inline function find_position_periodic( pos, k, boxsize )
    
    # re-define here to avoid type-instability
    box   = eltype(pos[1])(boxsize)
    
    x = (k & 0x1) == 0 ? pos[1] : add_subtr_box(pos[1], box)
	
	y = (k & 0x2) == 0 ? pos[2] : add_subtr_box(pos[2], box)

	z = (k & 0x4) == 0 ? pos[3] : add_subtr_box(pos[3], box)
            
    return x, y, z
end


@inline function find_position_periodic!( _pos::Array{T}, pos::Array{T}, k::Integer, boxsize::Float64) where T
    
    # re-define here to avoid type-instability
    box   = eltype(pos[1])(boxsize)

    _pos[1] = (k & 0x1) == 0 ? pos[1] : add_subtr_box(pos[1], box)
	
	_pos[2] = (k & 0x2) == 0 ? pos[2] : add_subtr_box(pos[2], box)

	_pos[3] = (k & 0x4) == 0 ? pos[3] : add_subtr_box(pos[3], box)

    return _pos
end