"""
irradiance.py
=============
Implements:
  - View factor surface→Sun   (Equation 30)
  - View factor surface→Earth (Equation 33, Richmond 2010 [Ref 57])
  - Albedo phase function psi (Equation 34)
  - Solar irradiance          (Equation 29)
  - Albedo irradiance         (Equation 32)
  - Earth IR irradiance       (Equation 35)
  - Total irradiance          (Equation 36)

All equations from:
    Morsch Filho et al. (2020). Energies 13(24), 6691.
"""

import numpy as np
from constants import RE, G_SUN, G_EARTH, ALBEDO


# =============================================================================
# VIEW FACTOR: SURFACE → SUN (Equation 30)
# =============================================================================

def view_factor_sun(Nk, r_sun_vec):
    """
    View factor of face k toward the Sun.
    Equation 30: F_k→sun = Nk · (r_sun / |r_sun|)

    Returns max(0, dot) — negative means face is away from Sun.
    """
    r_sun_hat = r_sun_vec / np.linalg.norm(r_sun_vec)
    return max(0.0, float(np.dot(Nk, r_sun_hat)))


# =============================================================================
# VIEW FACTOR: SURFACE → EARTH (Equation 33)
# Source: Richmond (2010) [Ref 57], analytical flat-plate/sphere formula
# =============================================================================

def view_factor_earth(Nk, r_vec):
    """
    View factor of face k toward the Earth disk.
    Equation 33a–f of paper, from Richmond (2010) [Ref 57].

    Parameters
    ----------
    Nk    : ndarray (3,) — face outward normal (unit vector)
    r_vec : ndarray (3,) — satellite position [km]
    """
    r = np.linalg.norm(r_vec)

    # Eq. 33d
    H = r / RE
    # Eq. 33c
    phi = np.arcsin(np.clip(1.0 / H, -1, 1))

    # Eq. 33b — note negative sign (inward component toward Earth)
    Nk_norm = np.linalg.norm(Nk)
    cos_gk = float(np.dot(-Nk, r_vec)) / (Nk_norm * r)
    cos_gk = np.clip(cos_gk, -1.0, 1.0)
    gamma_k = np.arccos(cos_gk)

    # Eq. 33a — three cases
    if gamma_k <= np.pi / 2.0 - phi:
        # Case 1: face fully sees Earth
        return np.cos(gamma_k) / H**2

    elif gamma_k <= np.pi / 2.0 + phi:
        # Case 2: partial view
        H2 = H**2
        sin_gk = np.sin(gamma_k)
        if abs(sin_gk) < 1e-12:
            return 0.0

        # Eq. 33e
        arg_W1 = np.clip(np.sqrt(H2 - 1.0) / (H * sin_gk), -1.0, 1.0)
        W1 = 0.5 * np.arcsin(arg_W1)

        # Eq. 33f
        tan_gk = sin_gk / cos_gk if abs(cos_gk) > 1e-12 else np.sign(sin_gk) * 1e12
        arg_acos = np.clip(-np.sqrt(H2 - 1.0) / tan_gk, -1.0, 1.0)
        inner = 1.0 - H2 * cos_gk**2
        if inner < 0:
            return 0.0
        W2 = (cos_gk * np.arccos(arg_acos)
              - np.sqrt(H2 - 1.0) * np.sqrt(inner)) / (2.0 * H2)

        return max(0.0, 0.5 - (2.0 / np.pi) * (W1 - W2))

    else:
        # Case 3: face points away from Earth
        return 0.0


# =============================================================================
# ALBEDO PHASE FUNCTION PSI (Equation 34)
# =============================================================================

def albedo_psi(r_vec, r_sun_vec):
    """
    Albedo reflection phase function.
    Equation 34: psi = b*cos(chi) for chi in [0, pi/2], else 0.

    chi is the angle between satellite position and Sun vectors
    (satellite's angular distance from subsolar point).
    """
    r_mag    = np.linalg.norm(r_vec)
    rsun_mag = np.linalg.norm(r_sun_vec)
    chi = np.arccos(np.clip(
        np.dot(r_vec, r_sun_vec) / (r_mag * rsun_mag), -1.0, 1.0))
    if chi <= np.pi / 2.0:
        return ALBEDO * np.cos(chi)   # Eq. 34
    return 0.0


# =============================================================================
# IRRADIANCE ON ALL FACES (Equations 29, 32, 35, 36)
# =============================================================================

def compute_irradiance(normals, r_vec, r_sun_vec, xi, A_faces, psi):
    """
    Compute irradiance [W] on each of the 6 faces and totals.

    Parameters
    ----------
    normals   : ndarray (6,3) — face normals in geocentric frame
    r_vec     : ndarray (3,) — satellite position [km]
    r_sun_vec : ndarray (3,) — Sun position [km]
    xi        : float — eclipse factor (0=eclipse, 1=sunlit)  Eq. 29
    A_faces   : list[float] — area of each face [m^2], length 6
    psi       : float — albedo phase function value

    Returns dict with per-face and total irradiance [W] and flux [W/m^2].
    """
    Q_sun  = np.zeros(6)
    Q_alb  = np.zeros(6)
    Q_e    = np.zeros(6)

    for k, (Nk, Ak) in enumerate(zip(normals, A_faces)):
        Fk_sun   = view_factor_sun(Nk, r_sun_vec)          # Eq. 30
        Fk_earth = view_factor_earth(Nk, r_vec)             # Eq. 33

        Q_sun[k] = G_SUN   * Ak * Fk_sun   * xi            # Eq. 29
        Q_alb[k] = G_SUN   * Ak * Fk_earth * psi * xi      # Eq. 32
        Q_e[k]   = G_EARTH * Ak * Fk_earth                  # Eq. 35

    Q_tot_W = Q_sun.sum() + Q_alb.sum() + Q_e.sum()         # Eq. 36
    A_total  = sum(A_faces)

    return {
        'Q_sun_W':  Q_sun.sum(),
        'Q_alb_W':  Q_alb.sum(),
        'Q_e_W':    Q_e.sum(),
        'Q_tot_W':  Q_tot_W,
        'Q_sun':    Q_sun.sum() / A_total,    # flux W/m^2
        'Q_alb':    Q_alb.sum() / A_total,
        'Q_e':      Q_e.sum()   / A_total,
        'Q_tot':    Q_tot_W     / A_total,
        'Q_faces':  Q_sun + Q_alb + Q_e,      # per-face total W
        'Q_faces_flux': (Q_sun + Q_alb + Q_e) / np.array(A_faces),  # W/m^2
    }


# =============================================================================
# ORBIT-AVERAGED IRRADIANCE FOR FIGURE 6
# =============================================================================

def orbit_avg_irradiance(beta_deg, alt_km):
    """
    Compute orbit-averaged irradiance flux [W/m^2] for Figure 6.
    Uses analytical eclipse fraction (Eqs. 40–41) and average view factors.

    The paper plots average flux incident on a unit area facing the source.
    For a sphere at distance r, the average view factor to Earth is:
        F_earth_avg = 0.5 * (1 - sqrt(1 - (RE/r)^2))
    This is the standard result for a flat plate parallel to local horizon
    averaged over a full orbit (Richmond 2010 [Ref 57]).
    """
    from orbital_mechanics import eclipse_fraction

    r = RE + alt_km
    fE = eclipse_fraction(beta_deg, alt_km)

    # Average Earth view factor (analytical, for isotropically oriented plate)
    F_earth_avg = 0.5 * (1.0 - np.sqrt(1.0 - (RE / r)**2))

    # Average solar flux over full orbit (sunlit fraction × G_SUN)
    Q_sun_avg = G_SUN * (1.0 - fE)

    # Albedo: average chi over sunlit arc. For beta=0, chi ranges 0→90°.
    # Average of b*cos(chi) over sunlit arc ≈ b * (2/pi) for beta→0.
    # More precisely, use average cos(chi) = (2/pi)*sin(pi/2*(1-fE))
    # This is derived from integrating cos(chi) over the sunlit portion.
    sunlit_frac = 1.0 - fE
    if sunlit_frac > 0:
        avg_psi = ALBEDO * (2.0 / np.pi) * np.sin(np.pi / 2.0 * sunlit_frac)
    else:
        avg_psi = 0.0
    Q_alb_avg = G_SUN * F_earth_avg * avg_psi

    # Earth IR: constant (no eclipse dependence for thermal emission)
    Q_e_avg = G_EARTH * F_earth_avg

    return Q_sun_avg, Q_alb_avg, Q_e_avg, fE
