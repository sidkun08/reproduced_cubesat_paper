"""
sun_shadow.py
=============
Implements:
  - Sun position vector  (Equation 21)
  - Earth shadow check   (Equation 31, cylindrical shadow)

All equations from:
    Morsch Filho et al. (2020). Energies 13(24), 6691.
"""

import numpy as np
from constants import RE, AU


# =============================================================================
# SUN POSITION (Equation 21)
# =============================================================================

def sun_position(year, month, day, UT=0.0):
    """
    Compute Sun position in geocentric equatorial frame.
    Equations 21a–i of paper, from Curtis (2014) [Ref 40].

    Parameters
    ----------
    year, month, day : int
    UT               : float — Universal Time [hours from midnight]

    Returns
    -------
    u_hat      : ndarray (3,) — unit vector Earth→Sun
    r_sun_vec  : ndarray (3,) — Sun position vector [km]
    r_sun_mag  : float        — Sun–Earth distance [km]
    lam        : float        — ecliptic longitude [rad]  (= Gamma in Eq. 39)
    eps        : float        — obliquity [rad]
    """
    y, m, d = year, month, day

    # Eq. 21a: Julian Day Number
    JD = (367 * y
          - int(7 * (y + int((m + 9) / 12)) / 4)
          + int(275 * m / 9)
          + d + UT / 24.0
          + 1721013.5)

    # Eq. 21b: days since J2000
    nd = JD - 2451545.0

    # Eq. 21c: mean anomaly of Sun
    Ms = np.deg2rad(357.529 + 0.98560023 * nd)

    # Eq. 21d: mean longitude
    Ls = np.deg2rad(280.459 + 0.98564736 * nd)

    # Eq. 21e: ecliptic longitude
    lam = Ls + np.deg2rad(1.915) * np.sin(Ms) + np.deg2rad(0.020) * np.sin(2 * Ms)

    # Eq. 21f: obliquity
    eps = np.deg2rad(23.439 - 3.56e-7 * nd)

    # Eq. 21g: unit vector
    u_hat = np.array([
        np.cos(lam),
        np.cos(eps) * np.sin(lam),
        np.sin(eps) * np.sin(lam),
    ])

    # Eq. 21h: Sun–Earth distance [AU → km]
    r_sun_mag = (1.00014 - 0.01671 * np.cos(Ms) - 0.000140 * np.cos(2 * Ms)) * AU

    # Eq. 21i
    r_sun_vec = r_sun_mag * u_hat

    return u_hat, r_sun_vec, r_sun_mag, lam, eps


# =============================================================================
# SHADOW CHECK (Equation 31, cylindrical shadow)
# =============================================================================

def check_shadow(r_vec, r_sun_vec):
    """
    Determine eclipse state using cylindrical shadow model.
    Equations 31a–d of paper.

    Parameters
    ----------
    r_vec     : ndarray (3,) — satellite position [km]
    r_sun_vec : ndarray (3,) — Sun position [km]

    Returns
    -------
    xi : float — 0.0 (eclipse) or 1.0 (sunlit)
    """
    r_mag    = np.linalg.norm(r_vec)
    rsun_mag = np.linalg.norm(r_sun_vec)

    # Eq. 31a
    chi_c   = np.arccos(np.clip(RE / r_mag,    -1, 1))
    # Eq. 31b
    chi_sun = np.arccos(np.clip(RE / rsun_mag, -1, 1))
    # Eq. 31c
    chi     = np.arccos(np.clip(
        np.dot(r_vec, r_sun_vec) / (r_mag * rsun_mag), -1, 1))

    # Eq. 31d
    return 0.0 if (chi_c + chi_sun) <= chi else 1.0
