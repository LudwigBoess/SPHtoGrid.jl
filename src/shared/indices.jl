"""
    function calculate_index(i::Integer, j::Integer, x_pixels::Integer)

Calculates the index of a flattened 2D image array.
"""
function calculate_index(i::T, j::T, x_pixels::T) where {T}
    return floor(T, i * x_pixels + j) + 1
end

"""
    function calculate_index(i::Integer, j::Integer, x_pixels::Integer)

Calculates the index of a flattened 3D image array.
"""
function calculate_index(i::T, j::T, k::T, x_pixels::T, y_pixels::T) where {T}
    return floor(T, i * x_pixels * y_pixels + j * y_pixels + k) + 1
end