"""
thermal.py
==========
Steady-state temperature model (Equations 37–38).

All equations from:
    Morsch Filho et al. (2020). Energies 13(24), 6691.
"""

import numpy as np
from constants import SIGMA


def steady_state_temperature(Q_tot_W, A_total_m2, alpha_over_e=1.0):
    """
    Steady-state single-node temperature [K].

    Equation 38:  T = (alpha/e * Q_tot / (sigma * A))^(1/4)

    Derived from Eq. 37 thermal balance:
        alpha * Q_tot = e * sigma * A * (T^4 - T_inf^4)
    T_inf = 2.7 K (outer space) is negligible (stated in paper after Eq. 37).

    Parameters
    ----------
    Q_tot_W     : float — total absorbed irradiance [W]
    A_total_m2  : float — total external surface area [m^2]
    alpha_over_e: float — absorptivity / emissivity ratio

    Returns
    -------
    T : float — temperature [K]
    """
    if Q_tot_W <= 0.0:
        return 0.0
    T = (alpha_over_e * Q_tot_W / (SIGMA * A_total_m2)) ** 0.25
    return T


def orbit_avg_temperature(Q_sun_avg, Q_alb_avg, Q_e_avg,
                          A_total_m2, alpha_over_e, half_area=False):
    """
    Average temperature from orbit-averaged irradiance fluxes [W/m^2].

    For Figure 7: paper assumes half the external surface absorbs radiation
    (stated: "hypothetical geometry whose half of the external surface could
    absorb the incoming radiation"). So absorbed power = flux * A_total/2.

    For Figure 12: full attitude-dependent view factors used, not this function.

    Parameters
    ----------
    Q_sun_avg, Q_alb_avg, Q_e_avg : float — average flux [W/m^2]
    A_total_m2  : float — total external surface area [m^2]
    alpha_over_e: float
    half_area   : bool — if True, use A_total/2 as absorbing area (for Fig 7)
    """
    A_eff = A_total_m2 / 2.0 if half_area else A_total_m2
    Q_flux_total = Q_sun_avg + Q_alb_avg + Q_e_avg   # W/m^2
    Q_tot_W = Q_flux_total * A_eff                    # W
    return steady_state_temperature(Q_tot_W, A_total_m2, alpha_over_e)
