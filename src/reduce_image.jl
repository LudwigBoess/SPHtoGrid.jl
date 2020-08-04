
"""
    function reduce_image_2D( image::Array{<:Real}, w_image::Array{<:Real},
                                            x_pixels::Int64, y_pixels::Int64)

Unflattens an image array to a 2D array.
"""
@inline @fastmath function reduce_image_2D( image::Array{<:Real}, w_image::Array{<:Real},
                                            x_pixels::Int64, y_pixels::Int64)


    im_plot = zeros(x_pixels, y_pixels)
    k = 1
    @inbounds for i = 1:x_pixels, j = 1:y_pixels
        im_plot[i,j] = image[k]
        if w_image[k] > 0.0
            im_plot[i, j] /= w_image[k]
        end
        k += 1
    end
    return im_plot
end


@inline @fastmath function reduce_image_2D( image::Array{<:Real},
                                            x_pixels::Int64, y_pixels::Int64)


    im_plot = zeros(x_pixels, y_pixels)
    k = 1
    @inbounds for i = 1:x_pixels, j = 1:y_pixels
        im_plot[i,j] = image[k,1]
        if image[k,2] > 0.0
            im_plot[i, j] /= image[k,2]
        end
        k += 1
    end
    return im_plot
end

"""
    function reduce_image_3D( image::Array{<:Real}, w_image::Array{<:Real},
                                            x_pixels::Int64, y_pixels::Int64, z_pixels::Int64)

Unflattens an image array to a 3D array.
"""
@inline @fastmath function reduce_image_3D( image::Array{<:Real},
                                            x_pixels::Int64, y_pixels::Int64, z_pixels::Int64)


    im_plot = zeros(x_pixels, y_pixels, z_pixels)
    m = 1
    @inbounds for i = 1:x_pixels, j = 1:y_pixels, k = 1:z_pixels
        im_plot[i,j,k] = image[m,1]
        if image[m,1] > 0.0
            im_plot[i, j, k] /= image[m,2]
        end
        m += 1
    end
    return im_plot
end

@inline @fastmath function reduce_image_3D( image::Array{<:Real}, w_image::Array{<:Real},
                                         x_pixels::Int64, y_pixels::Int64, z_pixels::Int64)


    im_plot = zeros(x_pixels, y_pixels, z_pixels)
    m = 1
    @inbounds for i = 1:x_pixels, j = 1:y_pixels, k = 1:z_pixels
        im_plot[i,j,k] = image[m]
        if w_image[m] > 0.0
            im_plot[i, j, k] /= w_image[m]
        end
        m += 1
    end
    return im_plot
end

"""
    function reduce_futures(fut::Array{<:Tuple})

Reduces the touple returned by the Array of `Future`s to image arrays. 
"""
@inline @fastmath function reduce_futures(fut::Array{<:Tuple})
    
    image   = fut[1][1]
    w_image = fut[1][2]

    N = length(fut)
    if N > 1
        @inbounds for i = 1:N
            image   .+= fut[i][1]
            w_image .+= fut[i][2]
        end
    end
    return image, w_image
end