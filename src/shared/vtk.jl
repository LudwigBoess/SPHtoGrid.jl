using WriteVTK

"""
    write_vtk_image(filename::String, image::Array{<:Real},
                    image_name::String, par::mappingParameters;
                    units::String = "[i.u.]", snap::Integer = 0)

Writes a mapped image to a vti file.
"""
function write_vtk_image(filename::String, image::Array{<:Real},
                        image_name::String, par::mappingParameters;
                        units::String = "[i.u.]", snap::Integer = 0)

    # construct pixel centers
    x, y, z = get_map_grid_3D(par)

    vtk_grid(filename, x, y, z) do vtk
        vtk[image_name] = image  # scalar field attached to points
        vtk["Units"] = units     # metadata ("field data" in VTK)
        vtk["Snap"] = snap       # metadata ("field data" in VTK)
    end
end