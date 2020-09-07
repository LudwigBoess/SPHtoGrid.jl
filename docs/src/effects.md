# Effects

If you want to map more than the primitive output quantities to a grid you can use some helper functions to map special effects.

## 2D surface density

To get the projected 2D surface density from your data you can use

```julia
density_2D(rho::Real, pixelSideLength::Real, 
           Mass::Real=1.989e43, Length::Real=3.085678e21)
```

where `rho` is the density in internal units, `pixelSideLength` can be taken from the `mappingParameters` `struct` and `Mass` and `Length` define the mass and length unit of the simulation respectively.
To use this function you need to map the quantity without weighting, so:

```julia
par     = mappingParameters(...)
bin_q   = density_2D.(rho, par.pixelSideLength)
weights = part_weight_one(length(bin_q))
```

## X-Ray emission

To compute the contribution of each particle to the total (thermal) X-ray emission you can use

```julia
x_ray_emission( n_cm3::Real, T_eV::Real; 
                Emin::Real=5.e4, Emax::Real=1.e10, 
                xH::Real=0.76)
```
where `n_cm3` is the particle density in ``1/cm^3``, `T_eV` is the temperature in electron volts. The optional input arguments `Emin` and `Emax` give the minimum and maximum energy range of the observational instrument, while `xH` gives the hydrogen fraction of the simulation.

The correct weight function in this case is

```julia
weights = part_weight_physical(length(bin_q), par)
```

## Sunyaev-Z'eldovich Effect

### Thermal

To compute the contribution of each particle to the thermal SZ-effect you can use

```julia
thermal_SZ( n_cm3::Real, T::Real, 
            z::Real=0.0, ν::Real=1.44e9; 
            DI_over_I::Bool=false )
```

where `n_cm3` is the particle density in ``1/cm^3``, `T` is the temperature in Kelvin, `z` is the redshift and `ν` is the observation frequency. `DI_over_I` outputs in units of ``dI/I`` if set to `true` and `dT/T` otherwise.

### Kinetic

To compute the contribution of each particle to the kinetic SZ-effect you can use

kinetic_SZ(n_cm3::Real, vel_y_cgs::Real, 
                    ν::Real=1.e9, z::Real=0.0; 
                    DI_over_I::Bool=false)

where `n_cm3` is the particle density in ``1/cm^3``, `vel_y_cgs` is the velocity in y-direction to the projection plane in `cm/s`, `z` is the redshift and `ν` is the observation frequency. `DI_over_I` outputs in units of ``dI/I`` if set to `true` and `dT/T` otherwise.

### Weights
Independent of what version of the SZ effect you want to map you need to use the physical weight function

```julia
weights = part_weight_physical(length(bin_q), par)
```