"""
orbital_mechanics.py
====================
Implements:
  - TLE → COE  (Equations 1–4)
  - COE → position  (Equations 5–8)
  - Drag perturbation  (Equations 12–17)
  - J2 perturbation  (Equation 18)
  - Gauss variational equations  (Equation 19)
  - RK4 integrator
  - Beta angle  (Equation 39)
  - Eclipse fraction  (Equations 40–41)

All equations from:
    Morsch Filho et al. (2020). Energies 13(24), 6691.
    https://doi.org/10.3390/en13246691
"""

import numpy as np
from constants import MU, RE, J2, CD


# =============================================================================
# ROTATION MATRICES (Equations 7a, 7b, 7c)
# =============================================================================

def R1(a):
    """Rotation matrix about x-axis (Eq. 7b)."""
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0],
                     [0, c, s],
                     [0,-s, c]], dtype=float)

def R2(a):
    """Rotation matrix about y-axis."""
    c, s = np.cos(a), np.sin(a)
    return np.array([[ c, 0,-s],
                     [ 0, 1, 0],
                     [ s, 0, c]], dtype=float)

def R3(a):
    """Rotation matrix about z-axis (Eqs. 7a, 7c)."""
    c, s = np.cos(a), np.sin(a)
    return np.array([[ c, s, 0],
                     [-s, c, 0],
                     [ 0, 0, 1]], dtype=float)


# =============================================================================
# KEPLER EQUATION SOLVER (Equation 3)
# =============================================================================

def solve_kepler(M, e, tol=1e-12, max_iter=200):
    """
    Solve E - e*sin(E) = M for eccentric anomaly E using Newton-Raphson.
    Eq. 3: E - e*sin(E) = M
    """
    E = float(M)
    for _ in range(max_iter):
        dE = (M - E + e * np.sin(E)) / (1.0 - e * np.cos(E))
        E += dE
        if abs(dE) < tol:
            break
    return E


# =============================================================================
# TLE → CLASSICAL ORBITAL ELEMENTS (Equations 1–4)
# =============================================================================

def tle_to_coe(i_deg, Omega_deg, e, omega_deg, M_deg, n_rev_day):
    """
    Convert TLE parameters to Classical Orbital Elements.

    Equations 1–4 of the paper.

    Returns dict: {h [km^2/s], e [-], Omega [rad], i [rad], omega [rad], theta [rad]}
    """
    # Mean motion: convert rev/day → rad/s
    n = n_rev_day * 2.0 * np.pi / 86400.0          # rad/s

    # Eq. 1: semi-major axis
    a = (MU / n**2) ** (1.0 / 3.0)                 # km

    # Eq. 2: specific angular momentum
    h = np.sqrt(a * MU * (1.0 - e**2))             # km^2/s

    # Eq. 3: eccentric anomaly (Newton-Raphson)
    M_rad = np.deg2rad(M_deg)
    E = solve_kepler(M_rad, e)

    # Eq. 4: true anomaly
    theta = 2.0 * np.arctan2(
        np.sqrt(1.0 + e) * np.sin(E / 2.0),
        np.sqrt(1.0 - e) * np.cos(E / 2.0)
    )

    return {
        'h':     h,
        'e':     e,
        'Omega': np.deg2rad(Omega_deg),
        'i':     np.deg2rad(i_deg),
        'omega': np.deg2rad(omega_deg),
        'theta': theta,
    }


# =============================================================================
# COE → POSITION VECTOR (Equations 5–8)
# =============================================================================

def coe_to_position(coe):
    """
    Convert COE to position vector in geocentric equatorial frame.

    Equations 5–8 of the paper.

    Returns
    -------
    r_vec : ndarray (3,) — position [km]
    r_mag : float        — orbital radius [km]
    Q     : ndarray (3,3) — rotation matrix geocentric→perifocal
    """
    h, e, Omega, i, omega, theta = (
        coe['h'], coe['e'], coe['Omega'],
        coe['i'], coe['omega'], coe['theta']
    )

    # Eq. 5: position in perifocal frame
    r_mag = h**2 / MU / (1.0 + e * np.cos(theta))
    r_peri = r_mag * np.array([np.cos(theta), np.sin(theta), 0.0])

    # Eq. 6–7: rotation matrix Q = R3(omega) @ R1(i) @ R3(Omega)
    Q = R3(omega) @ R1(i) @ R3(Omega)

    # Eq. 8: r = Q^T @ r_peri  (Q orthogonal → inverse = transpose)
    r_vec = Q.T @ r_peri

    return r_vec, r_mag, Q


def coe_radius(coe):
    """Orbital radius [km] from COE."""
    return coe['h']**2 / MU / (1.0 + coe['e'] * np.cos(coe['theta']))


def orbital_velocity_mag(coe):
    """Orbital speed magnitude [km/s] — Eq. 16."""
    h, e, theta = coe['h'], coe['e'], coe['theta']
    return MU / h * np.sqrt(np.sin(theta)**2 + (e + np.cos(theta))**2)


def coe_apogee_perigee(coe):
    """Return (apogee_alt_km, perigee_alt_km)."""
    h, e = coe['h'], coe['e']
    a = h**2 / MU / (1.0 - e**2)
    r_a = a * (1 + e)
    r_p = a * (1 - e)
    return r_a - RE, r_p - RE


# =============================================================================
# DRAG PERTURBATION (Equations 12–17)
# =============================================================================

def drag_perturbation(coe, rho_kg_m3, cd=CD, A_over_m=0.01654):
    """
    Atmospheric drag perturbation components in LVLH frame [km/s^2].

    Parameters
    ----------
    coe        : dict of COE
    rho_kg_m3  : float — atmospheric density [kg/m^3]
    cd         : float — drag coefficient (default 2.2, Section 3.1.1)
    A_over_m   : float — frontal area / mass [m^2/kg]
                 Default: CD*A/m = 2.2*0.01/1.33 for 1U CubeSat

    Returns (prD, psD, pwD) in km/s^2.

    Equations 12–17 of paper.
    """
    h, e, theta = coe['h'], coe['e'], coe['theta']

    # Eq. 16: velocity magnitude [km/s]
    v_km = orbital_velocity_mag(coe)
    v_ms = v_km * 1e3                               # [m/s]

    # Eq. 15: drag scalar acceleration [m/s^2]
    pD_ms = -0.5 * rho_kg_m3 * v_ms**2 * (cd * A_over_m)
    pD_km = pD_ms * 1e-3                            # [km/s^2]

    # Orbital radius
    r = h**2 / MU / (1.0 + e * np.cos(theta))

    # Eq. 17a
    prD = (pD_km / v_km) * (MU * e * np.sin(theta) / h)
    # Eq. 17b
    psD = (pD_km / v_km) * (h / r)
    # Eq. 17c
    pwD = 0.0

    return prD, psD, pwD


# =============================================================================
# J2 GRAVITATIONAL PERTURBATION (Equation 18)
# =============================================================================

def j2_perturbation(coe):
    """
    J2 gravitational perturbation components in LVLH frame [km/s^2].
    Equation 18a–c of paper.
    """
    h, e, theta, i, omega = (
        coe['h'], coe['e'], coe['theta'], coe['i'], coe['omega']
    )
    r = h**2 / MU / (1.0 + e * np.cos(theta))
    c = -1.5 * J2 * MU * RE**2 / r**4

    prJ2 = c * (1.0 - 3.0 * np.sin(i)**2 * np.sin(omega + theta)**2)  # Eq. 18a
    psJ2 = c * np.sin(i)**2 * np.sin(2.0 * (omega + theta))            # Eq. 18b
    pwJ2 = c * np.sin(2.0 * i) * np.sin(omega + theta)                 # Eq. 18c

    return prJ2, psJ2, pwJ2


# =============================================================================
# GAUSS VARIATIONAL EQUATIONS (Equation 19)
# =============================================================================

def gauss_variational(coe, pr, ps, pw):
    """
    Time derivatives of COE under perturbation.
    Equations 19a–f of paper.

    Note: singularities at e=0 or i=0 (stated after Eq.19).
    PropCube-2 has e=0.00026, i=51.63° — no issue.
    """
    h, e, theta, Omega, i, omega = (
        coe['h'], coe['e'], coe['theta'],
        coe['Omega'], coe['i'], coe['omega']
    )
    r = h**2 / MU / (1.0 + e * np.cos(theta))

    dh     = r * ps                                                         # Eq. 19a
    de     = (h / MU * np.sin(theta) * pr
              + ps / (MU * h) * ((h**2 + MU * r) * np.cos(theta)
                                  + MU * e * r))                            # Eq. 19b
    dtheta = (h / r**2
              + h**2 * np.cos(theta) / (MU * e * h) * pr
              - (r + h**2 / MU) * np.sin(theta) / (e * h) * ps)            # Eq. 19c
    dOmega = r / (h * np.sin(i)) * np.sin(omega + theta) * pw              # Eq. 19d
    di     = r / h * np.cos(omega + theta) * pw                            # Eq. 19e
    domega = (-(h**2 * np.cos(theta)) / (MU * e * h) * pr
              + (r + h**2 / MU) * np.sin(theta) / (e * h) * ps
              - r * np.sin(omega + theta) / (h * np.tan(i)) * pw)          # Eq. 19f

    return {'dh': dh, 'de': de, 'dtheta': dtheta,
            'dOmega': dOmega, 'di': di, 'domega': domega}


def _add_coe(c, d, s):
    """Helper: c + s*d for RK4."""
    return {
        'h':     c['h']     + s * d['dh'],
        'e':     max(c['e'] + s * d['de'], 1e-7),
        'theta': c['theta'] + s * d['dtheta'],
        'Omega': c['Omega'] + s * d['dOmega'],
        'i':     c['i']     + s * d['di'],
        'omega': c['omega'] + s * d['domega'],
    }


def rk4_step(coe, dt, pr, ps, pw):
    """
    4th-order Runge-Kutta integration of Gauss variational equations.
    Perturbation components held constant over the step.
    """
    k1 = gauss_variational(coe, pr, ps, pw)
    k2 = gauss_variational(_add_coe(coe, k1, 0.5 * dt), pr, ps, pw)
    k3 = gauss_variational(_add_coe(coe, k2, 0.5 * dt), pr, ps, pw)
    k4 = gauss_variational(_add_coe(coe, k3, dt), pr, ps, pw)

    return {
        'h':     coe['h']     + dt / 6.0 * (k1['dh']     + 2*k2['dh']     + 2*k3['dh']     + k4['dh']),
        'e':     max(coe['e'] + dt / 6.0 * (k1['de']     + 2*k2['de']     + 2*k3['de']     + k4['de']), 1e-7),
        'theta': coe['theta'] + dt / 6.0 * (k1['dtheta'] + 2*k2['dtheta'] + 2*k3['dtheta'] + k4['dtheta']),
        'Omega': coe['Omega'] + dt / 6.0 * (k1['dOmega'] + 2*k2['dOmega'] + 2*k3['dOmega'] + k4['dOmega']),
        'i':     coe['i']     + dt / 6.0 * (k1['di']     + 2*k2['di']     + 2*k3['di']     + k4['di']),
        'omega': coe['omega'] + dt / 6.0 * (k1['domega'] + 2*k2['domega'] + 2*k3['domega'] + k4['domega']),
    }


# =============================================================================
# BETA ANGLE (Equation 39) and ECLIPSE FRACTION (Equations 40–41)
# =============================================================================

def beta_angle(coe, ecliptic_longitude_rad, obliquity_rad):
    """
    Orbit beta angle [deg] — Equation 39.

    Parameters
    ----------
    coe                   : dict of current COE
    ecliptic_longitude_rad: float — lambda from sun_position (Eq. 21e)
    obliquity_rad         : float — eps from sun_position (Eq. 21f)
    """
    Omega, i = coe['Omega'], coe['i']
    lam = ecliptic_longitude_rad
    eps = obliquity_rad

    beta = np.arcsin(
        np.cos(lam) * np.sin(Omega) * np.sin(i)
        - np.sin(lam) * np.cos(eps) * np.cos(Omega) * np.sin(i)
        + np.sin(lam) * np.sin(eps) * np.cos(i)
    )
    return np.rad2deg(beta)


def eclipse_fraction(beta_deg, alt_km):
    """
    Analytical eclipse fraction — Equations 40–41.

    Parameters
    ----------
    beta_deg : float — beta angle [degrees]
    alt_km   : float — orbital altitude [km]
    """
    r = RE + alt_km
    # Eq. 41
    beta_star_deg = np.rad2deg(np.arcsin(RE / r))

    if abs(beta_deg) >= beta_star_deg:
        return 0.0  # no eclipse

    # Eq. 40
    beta_r = np.deg2rad(beta_deg)
    numer = np.sqrt((r - RE)**2 + 2 * RE * (r - RE))
    denom = r * np.cos(beta_r)
    if denom <= 0:
        return 0.0
    fE = (1.0 / 180.0) * np.rad2deg(np.arccos(np.clip(numer / denom, -1, 1)))
    return np.clip(fE, 0.0, 1.0)
