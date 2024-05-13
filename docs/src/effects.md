# Examples

You can find some examples here for common maps. These examples assume that you also have `GadgetIO.jl`, `GadgetUnits.jl` and `SPHKernels.jl` installed.

In all examples `snap_base` points to the snapshot you want to map and `map_path` points to the folder where you want to store the fits files.

## Read data and convert units

First we need to read the data and convert it to physical units (unless you want to bother with comoving coordinates and `h0` in cosmological simulations).

Here `halo_pos` is the position of the halo you want to map and `rvir` is its virial radius.

```julia
# read the header
h = read_header(snap_base)
# define a unit conversion struct from `GadgetUnits.jl`
GU = GadgetPhysical(h)

# define the blocks you want to read
blocks = ["POS", "VEL", "HSML", "RHO", "U", "MASS", "BFLD", "MACH"]

# read all particles in a cubic volume around the halo
data = read_particles_in_volume(snap_base, blocks, halo_pos, rvir)

# convert to physical code units for mapping
pos  = data["POS"]  .* GU.x_physical
hsml = data["HSML"] .* GU.x_physical
rho  = data["RHO"]  .* GU.rho_physical
mass = data["MASS"] .* GU.m_physical

# we want to map the entire cube we read
xy_size = 2rvir
z_size  = 2rvir

# define mapping parameters and convert to physical units
par = mappingParameters(center = halo_pos .* GU.x_physical,
                        x_size = xy_size * GU.x_physical,
                        y_size = xy_size * GU.x_physical,
                        z_size = z_size * GU.x_physical,
                        Npixels = 1024,
                        boxsize=h.boxsize)

# define the kernel you want to use for mapping
k = WendlandC4(Float64, 2)
```

## 2D surface density

To get the projected 2D surface density from your data you can use

```julia
# you need to make a copy of the positions if you plan to re-use them
# the mapping shifts them in place
pos_map = copy(pos)

# convert density to physical cgs units
rho_gcm3 = data["RHO"] .* GU.rho_cgs

# to get the integrated values along the LOS you need physical weights and not reduce the image
weights = part_weight_physical(length(hsml), par, GU.x_cgs)

# actual mapping
map = sphMapping(pos_map, hsml, mass, rho,
                rho_gcm3, weights,
                param = par, kernel = k,
                reduce_image = false)

# filename of the output image
fo_image = map_path * "rho.fits"

# store the fits image
write_fits_image(fo_image, quantitiy_map, par, snap = snap, units = "g/cm^2")
```

The units in this case are of course ``g/cm^2``.


## Magnetic Field

To map the mean magnetic field along the LOS you need to use the density as weight (this is the default behaviour, so you can leave the arguemnt empty) and set `reduce_image=true`

```julia
# you need to make a copy of the positions if you plan to re-use them
# the mapping shifts them in place
pos_map = copy(pos)

# compute absolute value of magnetic field in muG
B = @. 1.e6 * √(data["BFLD"][1, :]^2 + data["BFLD"][2, :]^2 + data["BFLD"][3, :]^2)


# actual mapping
map = sphMapping(pos_map, hsml, mass, rho, B,
                param = par, kernel = k,
                reduce_image = true)

# filename of the output image
fo_image = map_path * "B.fits"

# store the fits image
write_fits_image(fo_image, quantitiy_map, par, snap = snap, units = "muG")
```

## Faraday Rotation

You can map the Faraday Rotation due to the LOS magnetic field strength for RM maps, or for continuous rotation of polarized synchrotron emission using
```@docs
rotation_measure
```
For the latter we provide the simple helper function
```@docs
get_dz
```

## X-Ray emission

To map the Xray surface brightness please use

```@docs
x_ray_emissivity
```

Here is an example code
```julia
# you need to make a copy of the positions if you plan to re-use them
# the mapping shifts them in place
pos_map = copy(pos)

# get temperature in keV
T_keV = get_T_keV(data["U"], data["MASS"], GU.T_eV)

# convert density to physical cgs units
rho_gcm3 = data["RHO"] .* GU.rho_cgs

# calculate X-ray emission per particle in the energy band Emin = 0.1 keV, Emax = 2.4 keV
Xray = x_ray_emissivity(T_keV, rho_cgs,
                        E0=0.1, E1=2.4)

# to get the integrated values along the LOS you need physical weights and not reduce the image
weights = part_weight_physical(length(hsml), par, GU.x_cgs)

# actual mapping
map = sphMapping(pos_map, hsml, mass, rho, Xray, weights,
                param = par, kernel = k,
                reduce_image = false)

# filename of the output image
fo_image = map_path * "Xray.fits"

# store the fits image
write_fits_image(fo_image, quantitiy_map, par, snap = snap, units = "erg/s/cm^2")
```

This returns a map in the units ``erg/s/cm^2``.

## Sunyaev-Z'eldovich Effect


You can compute compton-Y parameter, thermal and kinetic SZ effect with these functions:

```@docs
comptonY
thermal_SZ
kinetic_SZ
```

Example for the thermal SZ effect
```julia
# you need to make a copy of the positions if you plan to re-use them
# the mapping shifts them in place
pos_map = copy(pos)

# get temperature in K
T_K = data["U"] .* GU.T_K

# convert code density to electron density
n_cm3 = data["RHO"] .* GU.rho_ncm3

# calculate thermal SZ for redshift given in header
th_SZ = thermal_SZ(n_cm3, T_K, h.z)

# to get the integrated values along the LOS you need physical weights and not reduce the image
weights = part_weight_physical(length(hsml), par, GU.x_cgs)

# actual mapping
map = sphMapping(pos_map, hsml, mass, rho, th_SZ, weights,
                param = par, kernel = k,
                reduce_image = false)

# filename of the output image
fo_image = map_path * "ThermalSZ.fits"

# store the fits image
write_fits_image(fo_image, quantitiy_map, par, snap = snap, units = "")
```


## Gamma Ray emission

You can compute gamma-ray related maps as in [Pfrommer & Enßlin (2004)](https://ui.adsabs.harvard.edu/abs/2004A%26A...413...17P/abstract) with:

```@docs
jγ_PE04
λγ_PE04
gamma_luminosity_pions_PE04
gamma_flux_pions_PE04
```

Example for gamma-ray luminosity

```julia
# you need to make a copy of the positions if you plan to re-use them
# the mapping shifts them in place
pos_map = copy(pos)

# convert density to physical cgs units
rho_gcm3 = data["RHO"] .* GU.rho_cgs

# convert mass to physical cgs units
m_cgs = data["MASS"] .* GU.m_cgs

# get temperature in keV
T_keV = data["U"] .* GU.T_K

# calculate Iγ-ray luminosity per particle in the energy band Emin = 0.1 GeV, Emax = 200 GeV
# for proton spectra with slope 2.235
Lγ = gamma_luminosity_pions_PE04.(T_keV, m_cgs, rho_cgs, 2.235)

# weights ones means you sum up the values along the LOS
weights = ones(length(Lγ))

# actual mapping
map = sphMapping(pos_map, hsml, mass, rho, Lγ, weights,
                param = par, kernel = k,
                reduce_image = true)

# filename of the output image
fo_image = map_path * "L_gamma.fits"

# store the fits image
write_fits_image(fo_image, quantitiy_map, par, snap = snap, units = "GeV/s")
```


## Synchrotron Emission

I implemented a number of different functions to compute synchrotron emissivity by shock accelerated electrons.
Please use with caution as there are still some errors to be fixed!

```@docs
analytic_synchrotron
analytic_synchrotron_GS
analytic_synchrotron_Longair
analytic_synchrotron_HB07
```

## Stellar/DM density

If you want to map the mass density of a non-SPH particle, i.e. a particle distribution that does not have their density computed in the output snapshot you can do this with `mass_density`. This will perform an SPH loop to compute the mass density at particle position. Useful if you want to compute Stellar/DM density distributions.

Please note that this function uses a fixed number of neighbors (for simplicity) so the results will diverge from a proper SPH loop.

```@docs
mass_density
```


# Image Functions

## General

```@docs
surface_brightness_to_luminosity
```

## Synchrotron Specific

```@docs
synchrotron_SB_to_luminosity
total_synch_luminosity_from_SB
beam_in_kpc
convert_Pnu_map_to_mJy_beam
polarisation_fraction
polarisation_angle
```