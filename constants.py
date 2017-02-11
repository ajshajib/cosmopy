# speed of light
c_kms = 299792.458  # km/s
c = 29979245800  # cm/s

# Stefan-Boltzmann constant
sigma_sb = 5.670367e-5 # in cgs (erg/cm^2/s/K^4)

# Unit converseion
Mpc_to_cm = 3.085677581e24

# pi
PI = 3.14159265358979323846

# Gravitational constant
G = 6.67259e-8 # in cgs (cm^3 g^-1 s^-2)

# Radiation parameter over c^2 in cgs (g cm^-3 K^-4)
a_b_c2 = 4. * sigma_sb / c ** 3

# constant for critical density 3 / 8 pi G
critdens_const = 3. / (8. * PI * G)

# Omega_r for H_0 = 100 km/s/Mpc and T_CMB = 2.725 K, Omega_r = const * (T_CMB/2.725)^4 / h^2
omega_r_const = a_b_c2 / critdens_const / (1e7/Mpc_to_cm)**2 * 2.725**4