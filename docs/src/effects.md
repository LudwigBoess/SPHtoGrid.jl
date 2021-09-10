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
### Units

The units in this case are of course ``g/cm^2``.


## X-Ray emission

To compute the contribution of each particle to the total (thermal) X-ray emission you can use

```julia
x_ray_emission( n_cm3::Real, T_K::Real; 
                Emin::Real=5.e4, Emax::Real=1.e10, 
                xH::Real=0.76)
```
where `n_cm3` is the particle density in ``1/cm^3``, `T_K` is the temperature in Kelvin. The optional input arguments `Emin` and `Emax` give the minimum and maximum energy range of the observational instrument, while `xH` gives the hydrogen fraction of the simulation.

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

```julia
thermal_SZ( n_cm3::Real, T_K::Real, 
            z::Real=0.0, ν::Real=1.44e9; 
            DI_over_I::Bool=false )
```

where `n_cm3` is the particle density in ``1/cm^3``, `T` is the temperature in Kelvin, `z` is the redshift and `ν` is the observation frequency. `DI_over_I` outputs in units of ``dI/I`` if set to `true` and `dT/T` otherwise.

### Kinetic

To compute the contribution of each particle to the kinetic SZ-effect you can use

```julia
kinetic_SZ( n_cm3::Real, vel_y_cgs::Real, 
            ν::Real=1.e9, z::Real=0.0; 
            DI_over_I::Bool=false )
```

where `n_cm3` is the particle density in ``1/cm^3``, `vel_y_cgs` is the velocity in y-direction to the projection plane in `cm/s`, `z` is the redshift and `ν` is the observation frequency. `DI_over_I` outputs in units of ``dI/I`` if set to `true` and `dT/T` otherwise.

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

```julia
analytic_synchrotron_emission( rho_cgs::Array{<:Real}, B_cgs::Array{<:Real},
                               T_K::Array{<:Real}, Mach::Array{<:Real};
                               xH::Real=0.76, dsa_model::Integer=1, ν0::Real=1.44e9,
                               K_ep::Real=0.01,
                               integrate_pitch_angle::Bool=true,
                               convert_to_mJy::Bool=false )
```

where `rho_cgs` is the density in ``g/cm^3``, `B_cgs` is the magnetic field in Gauss, `T_K` is the temperature in Kelvin and `Mach` is the shock Mach number.
This returns an `Array` with the synchrotron emission per particle in units ``erg/cm^3/Hz/s``.
The keyword argument `xH` gives the hydrogen fraction of the simulation, if the simulation was run without a chemical model. This has to be extended to work with a chemical model as well.
`dsa_model` defines the Diffuse-Shock-Acceleration model which should be used. It takes the values `0...4`. These refer to [`KR07_acc`](@ref), [`KR13_acc`](@ref), [`Ryu19_acc`](@ref), [`CS14_acc`](@ref) and [`P16_acc`](@ref) respectively.
`ν0` is the reference frequency for the synchrotron emission. Please note that this needs to be blue-shifted in the case of cosmological simulations.
`K_ep` is the ratio between electron and proton acceleration.
You can integrate over all angles between velocity and magnetic field vector by setting `integrate_pitch_angle=true`.
If you want to get the result in ``mJy/cm`` instead of ``erg/cm^3/Hz/s``.

### Synchrotron Emission by explicitly integrating the Spectrum

You can also calculate the synchrotron emission following [Donnert et. al. (2016)](https://academic.oup.com/mnras/article/462/2/2014/2589941), Eq. 17.

```julia
spectral_synchrotron_emission(rho_cgs::Real, B_cgs::Real,
                              T_K::Real, Mach::Real;
                              xH::Real=0.76, dsa_model::Integer=1, 
                              ν0::Real=1.44e9,
                              K_ep::Real=0.01,
                              Emin::Real=5.0e4,
                              Emax::Real=1.e10,
                              p_inj::Real=0.1, # in [me*c]
                              integrate_pitch_angle::Bool=true,
                              convert_to_mJy::Bool=false,
                              N_sample_bins::Integer=128)

```

Similar to before `rho_cgs` is the density in ``g/cm^3``, `B_cgs` is the magnetic field in Gauss, `T_K` is the temperature in Kelvin and `Mach` is the shock Mach number.
This returns the synchrotron emission per particle in units ``erg/cm^3/Hz/s``, based on the number density of CRs in relation to the thermal gas. This relation is defined by the `dsa_model` and `K_ep`.

The keyword argument have the same meaning as in [`analytic_synchrotron_emission`](@ref) with the addition of:
`Emin` and `Emax` which give the energy range over which you wish to integrate.
`p_inj` defines the momentum at which the Maxwell-Boltzmann distrubution transitions to a power-law.
`N_sample_bins` gives the number of bins used for the integration of the spectrum.


### Synchrotron Emission from a pre-defined spectrum

If you want to calculate the synchrotron emission of a spectrum that you get from e.g. an external Fokker-Planck solver you can use [`spectral_synchrotron_emission`](@ref):

```julia
spectral_synchrotron_emission( n_p::Vector{<:Real}, 
                               p::Vector{<:Real},
                               B_cgs::Real;
                               ν0::Real=1.44e9,
                               integrate_pitch_angle::Bool=false,
                               convert_to_mJy::Bool=false )
```

where `n_p` is the number density in ``1/cm`` at the momenta `p` and `B_cgs` is again the absolute value of the magnetic field in Gauss.

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