"""
    density_2D( rho::Real, pixelSideLength::Real, 
                Mass::Real=1.989e43, Length::Real=3.085678e21)

Computes the density in units of [g/cm^2].

## Arguments:
- `rho`: SPH particle density in physical code units.
- `pixelSideLength`: length of pixel in physical units.
- `Mass`: Mass unit in [g].
- `Length`: Length unit in [cm]

## Mapping settings
- weight function: [`part_weight_one`](@ref)
- reduce image: `true`
"""
function density_2D(rho::Real, pixelSideLength::Real,
                    Mass::Real=1.989e43, Length::Real=3.085678e21)
    return rho * ( (Mass / Length^2 ) * pixelSideLength / Length )
end
