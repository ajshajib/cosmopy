import numpy as np
from scipy.integrate import quad
import constants as const


class Cosmology(object):
    def __init__(self, params={}):
        """

        :param params: dictionary containing cosmological params, if any value is not provided Planck15 values are used as defaults
        """
        if 'H_0' in params:
            self.H_0 = params['H_0'] #in km/s/Mpc
            self.h = self.H_0/100.
        else:
            self.H_0 = 67.31
            self.h = 0.6731

        if 'sum_m_nu' in params:
            self.sum_m_nu = params['sum_m_nu']
        else:
            self.sum_m_nu = 0.715  # in eV

        if 'N_eff' in params:
            self.N_eff = params['N_eff']
        else:
            self.N_eff = 3.13

        if 'w_0' in params:
            self.w_0 = params['w_0']
        else:
            self.w_0 = -1.

        if 'w_a' in params:
            self.w_a = params['w_a']
        else:
            self.w_a = 0.

        if 'omega_m' in params:
            self.omega_m = params['omega_m']
        else:
            self.omega_m = 0.315

        if 'omega_gamma' in params:
            self.omega_gamma = params['omega_gamma']
            self.T_cmb = (params['omega_gamma'] / const.omega_r_const * self.h**2)**(0.25) * 2.725
        else:
            if 'T_cmb' in params:
                self.T_cmb = params['T_cmb']  # CMB temperature in K
            else:
                self.T_cmb = 2.725  # CMB temperature in K
            self.omega_gamma =  const.omega_r_const  * (self.T_cmb/2.725)**4 / self.h**2

        self.T_nu = 0.7137658555036082 * self.T_cmb  # background neutrino temperature in K, the constant is (4/11)^(1/3)

        if 'omega_v' in params:
            self.omega_v = params['omega_v']
        else:
            self.omega_v = 1. - self.omega_m - self.omega_gamma * (1. + self.get_relative_nu_density(0.))

        if 'omega_k' in params:
            self.omega_k = params['omega_k']
        else:
            self.omega_k = 1. - self.omega_m - self.omega_v - self.omega_gamma * (1. + self.get_relative_nu_density(0.))

        if self.omega_m + self.omega_v + self.omega_k + self.omega_gamma * (1. + self.get_relative_nu_density(0.)) != 1.:
            print("Warning: the cosmological parameters are over-defined!")

    def print_cosmology(self):
        """

        :return: nothing, prints out the cosmological parameters
        """
        print("%.12f %.12f %.12f %.12f" % (self.omega_m, self.omega_gamma, self.omega_v, self.omega_k))

    def get_relative_nu_density(self, z):
        """

        :param z: redshift
        :return: relative neutrino density for given redshift z, Omega_r = Omega_gamma ( 1. + relative neutrino density )
        """
        y = 187. * (self.sum_m_nu / 94.e-3) / (1. + z)
        A = 0.3173
        p = 1.83

        prefactor = 0.22710731766  # 7/8 (4/11)^4/3

        return prefactor * self.N_eff * (1. + (A * y) ** p) ** (1. / p)

    def get_e(self, z):
        """

        :param z: redshift
        :return: E(z)
        """
        """# returns dyanmic w for DE
        def w(z):
            return w_0 + w_a * z / (1. + z)
        def func(x):
            return w(x) / (1. + x)"""

        if self.w_a == 0. and self.w_0 == -1.:
            return np.sqrt(self.omega_m * (1. + z)**3 + self.omega_k * (1. + z)**2 + self.omega_gamma * (1. + z)**4 * (
                1. + self.get_relative_nu_density(z)) + self.omega_v)
        else:
            return np.sqrt(self.omega_m * (1. + z)**3 + self.omega_k * (1. + z)**2 + self.omega_gamma * (1. + z)**4 * (
                1. + self.get_relative_nu_density(z)) + self.omega_v * (1 + z)**(3.*(1.+self.w_0+self.w_a)) * np.exp(-3.*self.w_a*z/(1.+z)))

    def get_comoving_distance(self, z):
        """

        :param z: redshift of the object
        :return: comoving distance to the object
        """
        def func(x):
            return 1/self.get_e(x)

        return const.c_kms * quad(func, 0., z)[0] / self.H_0

    def get_transverse_comoving_distance(self, z):
        """

        :param z: redshift of the object
        :return: transverse comoving distance to the ojbect, Hogg 1999
        """
        D_c = self.get_comoving_distance(z)
        D_H = const.c_kms/self.H_0

        if self.omega_k > 0:
            return D_H / np.sqrt(self.omega_k) * np.sinh(np.sqrt(self.omega_k) * D_c/D_H)
        elif self.omega_k == 0:
            return D_c
        else:
            return D_H / np.sqrt(-self.omega_k) * np.sin(np.sqrt(-self.omega_k) * D_c/D_H)

    def get_angular_diameter_distance(self, z, z0=0.):
        """

        :param z: redshift of the target
        :param z0: initial redshift, default is 0
        :return: angular diameter distance to between redshifts z0 and z
        """
        if z0 == 0.:
            return self.get_transverse_comoving_distance(z)/(1.+z)
        else:
            D_M1 = self.get_transverse_comoving_distance(z0)
            D_M2 = self.get_transverse_comoving_distance(z)
            D_H = const.c_kms/self.H_0
            return (D_M2 * np.sqrt( 1. + self.omega_k * D_M1**2 / D_H**2) - D_M1 * np.sqrt(1. + self.omega_k * D_M2**2 / D_H**2)) / (1.+z)