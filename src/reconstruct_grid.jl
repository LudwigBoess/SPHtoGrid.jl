"""
    get_map_grid_2D(par::mappingParameters)

Reconstruct the 2D grid used for mapping.
"""
function get_map_grid_2D(par::mappingParameters)

    x_grid = zeros(par.Npixels[1])
    y_grid = zeros(par.Npixels[2])

    for i = 1:par.Npixels[1]
        x_grid[i] = par.x_lim[1] + ( i - 0.5 ) * par.pixelSideLength
    end

    for i = 1:par.Npixels[2]
        y_grid[i] = par.y_lim[1] + ( i - 0.5 ) * par.pixelSideLength
    end

    return x_grid, y_grid
end

"""
    get_map_grid_3D(par::mappingParameters)

Reconstruct the 3D grid used for mapping.
"""
function get_map_grid_3D(par::mappingParameters)

    x_grid = zeros(par.Npixels[1])
    y_grid = zeros(par.Npixels[2])
    z_grid = zeros(par.Npixels[3])

    for i = 1:par.Npixels[1]
        x_grid[i] = par.x_lim[1] + ( i - 0.5 ) * par.pixelSideLength
    end

    for i = 1:par.Npixels[2]
        y_grid[i] = par.y_lim[1] + ( i - 0.5 ) * par.pixelSideLength
    end

    for i = 1:par.Npixels[3]
        z_grid[i] = par.y_lim[1] + ( i - 0.5 ) * par.pixelSideLength
    end

    return x_grid, y_grid, z_grid
end