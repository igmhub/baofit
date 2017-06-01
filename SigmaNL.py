import numpy as np
from astropy.cosmology import Planck15

def SigmaNL(z=2.4, cosmology=Planck15, sigma8=0.831):
    """
    Calculate sigmas that approximate the Lagrangian displacement field
    for separations ~100 Mpc/h, following Seo & Eisenstein 2007,
    https://arxiv.org/abs/astro-ph/0701079.  See also equations (12,13)
    of https://arxiv.org/abs/astro-ph/0604361 and the accompanying text
    which describes the numerical simulations used to obtain the value
    12.4 Mpc/h used below.
    """
    # Tabulate log growth rate out to z.
    zz = np.linspace(0, z, 200)
    gamma = 0.55 + 0.05 * (1 + cosmology.w(1))
    f = cosmology.Om(zz) ** gamma
    # Calculate the growth function G(z) with G(0) = 0.758 so that
    # G(z) = 1/(1+z) at high z.
    G = 0.758 * np.exp(-np.trapz(f / (1 + zz), zz))
    # Rescale Sigma0 = 12.4 Mpc/h for sigma8.  Units are Mpc/h.
    Sigma0 = 12.4 * (sigma8 / 0.9)
    # Calculate Sigma_perp, Sigma_par
    Sigma_perp = Sigma0 * G
    Sigma_par = Sigma_perp * (1 + f[-1])

    print('perp: {0:.3f} Mpc/h, par: {1:.3f} Mpc/h'
          .format(Sigma_perp, Sigma_par))

if __name__ == '__main__':
    SigmaNL()
