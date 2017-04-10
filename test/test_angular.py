from astropy.cosmology import wCDM
from cosmopy import Cosmology
import numpy as np

n_sample = 300 # number of test samples
n_dim = 4

center = np.array([72., .5, .5, -1.5])  # H_0, omega_m, omega_v, w
scale = np.array([8., .5, .5, 1.]) # width to uniformly distribute cosmological parameters along one direction
params = center + np.random.uniform(low=-1., high=1., size=(n_sample, n_dim)) * scale

def isclose(a, b, rel_tol=1e-06, abs_tol=0.0):
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

for item in params:
    # reference distances from astropy
    cos = wCDM(H0=item[0], Om0=item[1], Ode0=item[2], w0=item[3], Tcmb0=0., Neff=0., m_nu=0.)
    D_A_ref = cos.angular_diameter_distance(0.5).value

    # distance computed from cosmopy
    param = {'H_0':cos.H0.value, 'omega_m':cos.Om0, 'omega_v':cos.Ode0, 'omega_gamma':cos.Ogamma0, 'omega_k':cos.Ok0, 'w_0':cos.w0, 'sum_m_nu':0., 'N_eff':0.}
    cosmo = Cosmology(param)
    D_A = cosmo.get_angular_diameter_distance(0.5)

    assert isclose(D_A, D_A_ref)