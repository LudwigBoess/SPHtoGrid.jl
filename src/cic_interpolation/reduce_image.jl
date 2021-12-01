
"""
    function reduce_image_2D( image::Array{<:Real},
                              x_pixels::Int64, y_pixels::Int64)

Unflattens an image array to a 2D array of pixels.
"""
@inline function reduce_image_2D( image::Matrix{<:Real},
                                            x_pixels::Int64, y_pixels::Int64)


    im_plot = zeros(y_pixels, x_pixels)
    k = 1
    @inbounds for i = 1:y_pixels, j = 1:x_pixels
        im_plot[j,i] = image[k,1]
        if image[k,2] > 0.0
            im_plot[j, i] /= image[k,2]
        end
        k += 1
    end
    return im_plot
end

"""
    function reduce_image_3D( image::Array{<:Real}, w_image::Array{<:Real},
                                            x_pixels::Int64, y_pixels::Int64, z_pixels::Int64)

Unflattens an image array to a 3D array of pixels.
"""
@inline @fastmath function reduce_image_3D( image::Matrix{<:Real},
                                            x_pixels::Int64, y_pixels::Int64, z_pixels::Int64)


    im_plot = zeros(z_pixels, y_pixels, x_pixels)
    m = 1
    @inbounds for i = 1:z_pixels, j = 1:y_pixels, k = 1:x_pixels

        im_plot[k,j,i] = image[m,1]
        
        if image[m,1] > 0.0
            im_plot[k, j, i] /= image[m,2]
        end
        m += 1
    end
    return im_plot
end
