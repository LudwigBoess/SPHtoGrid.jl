"""
    rotate_to_xz_plane!(par::mappingParameters)

Rotates an array of 3D positions into the xz-plane.
"""
function rotate_to_xz_plane(par::mappingParameters) 

    cen = [par.center[1], par.center[3], par.center[2]]

    xlim = copy(par.x_lim)
    ylim = copy(par.z_lim)
    zlim = copy(par.y_lim)
    
    return mappingParameters(center=cen, 
                             x_lim=xlim, y_lim=ylim, z_lim=zlim,
                             Npixels=maximum(par.Npixels), 
                             boxsize=par.boxsize)
end


"""
    rotate_to_yz_plane(par::mappingParameters)

Rotates an array of 3D positions into the yz-plane.
"""
function rotate_to_yz_plane(par::mappingParameters) 

    cen = [par.center[2], par.center[3], par.center[1]]

    xlim = copy(par.y_lim)
    ylim = copy(par.z_lim)
    zlim = copy(par.x_lim)
    
    return mappingParameters(center=cen,
                             x_lim=xlim, y_lim=ylim, z_lim=zlim,
                             Npixels=maximum(par.Npixels), 
                             boxsize=par.boxsize)
end