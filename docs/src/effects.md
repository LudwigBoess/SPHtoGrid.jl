# Effects

If you want to map more than the primitive output quantities to a grid you can use some helper functions to map special effects.

## 2D surface density

To get the projected 2D surface density from your data you can use

```@docs
density_2D
```

where `rho` is the density in internal units, `pixelSideLength` can be taken from the `mappingParameters` `struct` and `Mass` and `Length` define the mass and length unit of the simulation respectively.
To use this function you need to map the quantity without weighting, so:

```julia
par     = mappingParameters(...)
bin_q   = density_2D.(rho, par.pixelSideLength)
weights = part_weight_one(length(bin_q))
```
### Units

The units in this case are of course ``g/cm^2``.


## X-Ray emission

To compute the contribution of each particle to the total (thermal) X-ray emission you can use

```@docs
x_ray_emission
```

### Weights

The correct weight function in this case is

```julia
weights = part_weight_physical(size(bin_q,1), par, unit_lenght_in_cm)
```

### Units

This returns a map in the units ``erg/cm^2/s/Hz``.

## Sunyaev-Z'eldovich Effect

### Thermal

To compute the contribution of each particle to the thermal SZ-effect you can use

```@docs
thermal_SZ
```

where `n_cm3` is the particle density in ``1/cm^3``, `T` is the temperature in Kelvin, `z` is the redshift and `ν` is the observation frequency. `DI_over_I` outputs in units of ``dI/I`` if set to `true` and `dT/T` otherwise.

### Kinetic

To compute the contribution of each particle to the kinetic SZ-effect you can use

```@docs
kinetic_SZ
```


### Weights
Independent of what version of the SZ effect you want to map you need to use the physical weight function

```julia
weights = part_weight_physical(size(bin_q,1), par, unit_lenght_in_cm)
```

### Units

This produces a unitless map.


## Synchrotron Emission

### Analytic Synchrotron Emission

You can calculate the analytic synchrotron emission of a particle as described in Loungair 8.128 with

```@docs
analytic_synchrotron_emission
```


### Synchrotron Emission by explicitly integrating the Spectrum

You can also calculate the synchrotron emission following [Donnert et. al. (2016)](https://academic.oup.com/mnras/article/462/2/2014/2589941), Eq. 17.

```@docs
spectral_synchrotron_emission(::Real, ::Real, ::Real, ::Real)
```



### Synchrotron Emission from a pre-defined spectrum

If you want to calculate the synchrotron emission of a spectrum that you get from e.g. an external Fokker-Planck solver you can use [`spectral_synchrotron_emission`](@ref):

```@docs
spectral_synchrotron_emission(::Vector{<:Real}, ::Vector{<:Real}, ::Real)
```

### Weights

The correct weight function for all these functions is

```julia
weights = part_weight_physical(size(bin_q,1), par, unit_lenght_in_cm)
```

### Units

This gives a map in the units ``erg/cm^2/Hz/s`` which is equal to ``10^{26} mJy``.
To compute the total synchrotron emission of the image use:
```julia
# total radio emission of the image
J_ν = (par.Npixels[1] * par.pixelSideLength)^2 / par.Npixels[1]^2 * sum(image)

# convert to mJy if not done in the mapping
J_ν *= 1.e26
```