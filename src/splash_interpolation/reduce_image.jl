
"""
    function reduce_image_2D( image::Array{<:Real},
                              x_pixels::Int64, y_pixels::Int64)

Unflattens an image array to a 2D array of pixels.
"""
@inline @fastmath function reduce_image_2D(image::Vector{<:Real},
    x_pixels::Int64, y_pixels::Int64)


    im_plot = zeros(y_pixels, x_pixels)
    k = 1
    @inbounds for i = 1:y_pixels, j = 1:x_pixels
        im_plot[j, i] = image[k]
        k += 1
    end
    return im_plot
end

"""
    function reduce_image_3D( image::Array{<:Real}, w_image::Array{<:Real},
                                            x_pixels::Int64, y_pixels::Int64, z_pixels::Int64)

Unflattens an image array to a 3D array of pixels.
"""
@inline @fastmath function reduce_image_3D( image::Vector{<:Real},
                                            x_pixels::Int64, y_pixels::Int64, z_pixels::Int64)


    im_plot = zeros(z_pixels, y_pixels, x_pixels)
    m = 1
    @inbounds for i = 1:z_pixels, j = 1:y_pixels, k = 1:x_pixels

        im_plot[k,j,i] = image[m]
        m += 1
    end
    return im_plot
end
