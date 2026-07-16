
"""
    function reduce_image_2D( image::Array{<:Real},
                              x_pixels::Int64, y_pixels::Int64,
                              reduce_image::Bool )

Unflattens the deposited image into an `(x_pixels, y_pixels, N_images)` array.
The flattened deposition uses `idx = ix*y_pixels + iy + 1` (x is the outer,
y the inner axis), so the result is indexed `[ix, iy]` and works for
non-square maps (`x_pixels != y_pixels`).
"""
function reduce_image_2D( image::Matrix{<:Real},
                                  x_pixels::Int64, y_pixels::Int64,
                                  reduce_image::Bool)

    N_images = size(image, 2) - 1
    im_plot  = zeros(x_pixels, y_pixels, N_images)

    @inbounds for ix = 0:x_pixels-1, iy = 0:y_pixels-1
        k = ix * y_pixels + iy + 1
        for Nimage = 1:N_images
            val = image[k, Nimage]
            if reduce_image && (image[k, end] > 0.0)
                val /= image[k, end]
            end
            im_plot[ix+1, iy+1, Nimage] = val
        end
    end

    return im_plot
end

"""
    function reduce_image_3D( image::Array{<:Real},
                              x_pixels::Int64, y_pixels::Int64, z_pixels::Int64)

Unflattens the deposited cube into an `(x_pixels, y_pixels, z_pixels)` array.
The flattened deposition uses `idx = ix*y_pixels*z_pixels + iy*z_pixels + iz + 1`,
so the result is indexed `[ix, iy, iz]` and works for non-cubic maps.
"""
@inline @fastmath function reduce_image_3D( image::Matrix{<:Real},
                                            x_pixels::Int64, y_pixels::Int64, z_pixels::Int64)

    im_plot = zeros(x_pixels, y_pixels, z_pixels)

    @inbounds for ix = 0:x_pixels-1, iy = 0:y_pixels-1, iz = 0:z_pixels-1
        m   = ix * y_pixels * z_pixels + iy * z_pixels + iz + 1
        val = image[m, 1]
        if image[m, 2] > 0.0
            val /= image[m, 2]
        end
        im_plot[ix+1, iy+1, iz+1] = val
    end

    return im_plot
end
