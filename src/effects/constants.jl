global const k_B        = 1.38066e-16
global const c_light    = 2.9979e10
global const σ_T        = 6.65245e-25
global const m_e        = 9.10953e-28
global const m_p        = 1.6726e-24
const global E_p0       = 938.272088e-3             # proton mass in [GeV]
const global E_π0       = 134.9769e-3               # π^0 rest-mass in [GeV]
global const h_planck   = 6.6261e-27
global const q_e        = 1.602176487e-20 * c_light
global const yPrefac    = σ_T * k_B / m_e / c_light^2
global const eV2cgs     = 1.60218e-12
global const cgs2eV     = 1.0/eV2cgs

# synchrotron
const global C_crit = 3q_e / (4π * m_e^3 * c_light^5) # Donnert+16, MNRAS 462, 2014–2032 (2016), Eg. 20 
const global j_nu_prefac = √(3) * q_e^3 / (m_e * c_light^2)

global const γ_th          = 5.0/3.0
const global γ_cr = 4.0 / 3.0
global const mJy_factor = 1.e26       # conversion factor from [erg/cm^3/Hz/s] to mJy/cm.
global const f_p_prefac = 4π / c_light
global const C_j        = 2.42e-24   # factor for X-ray emissivity (BS1996, eq.4)
global const gg         = 1.2  # Gaunt factor used for X-ray emission
const global faraday_prefac = q_e^3 / (2π * m_e^2 * c_light^4) # Dennison 1980 prefac = 2.6e-17 [cgs]
const global kpc = 3.085678e21

