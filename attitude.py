"""
attitude.py
===========
Implements all 5 attitude models (Equations 22–28):
  - Nadir   (Section 3.3.1, Eqs. 23–24)
  - RAM     (Section 3.3.2, Eqs. 25–26)
  - Sun1    (Section 3.3.3, Eqs. 27–28, Theta_s1=0, Theta_s3=0)
  - Sun2    (Section 3.3.3, Eqs. 27–28, Theta_s1=0, Theta_s3=45°)
  - Sun3    (Section 3.3.3, Eqs. 27–28, Theta_s1=35.2643°, Theta_s3=45°)

All equations from:
    Morsch Filho et al. (2020). Energies 13(24), 6691.
"""

import numpy as np
from constants import MU, RAM_SPIN_REV_PER_ORBIT
from orbital_mechanics import R1, R2, R3

# =============================================================================
# INITIAL FACE NORMALS IN PERIFOCAL FRAME (Equation 22)
# =============================================================================

N0 = np.array([
    [-1, 0, 0],   # N1 — Eq. 22a
    [ 1, 0, 0],   # N2 — Eq. 22b
    [ 0,-1, 0],   # N3 — Eq. 22c
    [ 0, 1, 0],   # N4 — Eq. 22d
    [ 0, 0,-1],   # N5 — Eq. 22e
    [ 0, 0, 1],   # N6 — Eq. 22f
], dtype=float)


# =============================================================================
# RODRIGUES ROTATION AROUND ARBITRARY AXIS
# Used in Sun-fixed attitude (Section 3.3.3)
# =============================================================================

def rot_arbitrary(axis, angle):
    """Rotation matrix about arbitrary unit axis by angle (Rodrigues formula)."""
    ax = np.asarray(axis, dtype=float)
    n = np.linalg.norm(ax)
    if n < 1e-12:
        return np.eye(3)
    ax = ax / n
    c, s = np.cos(angle), np.sin(angle)
    ux, uy, uz = ax
    return np.array([
        [c + ux*ux*(1-c),     ux*uy*(1-c) - uz*s,  ux*uz*(1-c) + uy*s],
        [uy*ux*(1-c) + uz*s,  c + uy*uy*(1-c),      uy*uz*(1-c) - ux*s],
        [uz*ux*(1-c) - uy*s,  uz*uy*(1-c) + ux*s,   c + uz*uz*(1-c)   ],
    ])


# =============================================================================
# NADIR ATTITUDE (Equations 23–24)
# =============================================================================

def attitude_nadir(coe, Q):
    """
    Nadir-pointing: one face always toward Earth center.
    Rotates around w-hat by true anomaly theta (Eq. 23).

    Parameters
    ----------
    coe : dict — current COE (needs 'theta')
    Q   : ndarray (3,3) — geocentric→perifocal matrix from coe_to_position

    Returns
    -------
    normals : ndarray (6,3) — face normals in geocentric equatorial frame
    """
    theta = coe['theta']
    # Eq. 23: rotation around w-hat (z in perifocal) by theta
    Q_nad = R3(theta)
    # Eq. 24: Nk = Q^T @ Q_nad^T @ n0k
    Qt = Q.T  # perifocal → geocentric
    return np.array([Qt @ Q_nad.T @ n for n in N0])


# =============================================================================
# RAM ATTITUDE (Equations 25–26)
# =============================================================================

def attitude_ram(coe, Q, t_sec):
    """
    RAM attitude: long axis aligned with velocity, spinning at 4 rev/orbit.
    Section 3.3.2, Equations 25–26.

    Parameters
    ----------
    coe   : dict — current COE
    Q     : ndarray (3,3) — geocentric→perifocal matrix
    t_sec : float — elapsed simulation time [s]
    """
    h, e = coe['h'], coe['e']
    theta = coe['theta']

    # Orbital period
    a = h**2 / MU / (1.0 - e**2)
    T_orb = 2.0 * np.pi * np.sqrt(a**3 / MU)   # seconds

    # Spin angle at time t (4 rev/orbit per paper Section 4)
    Theta = RAM_SPIN_REV_PER_ORBIT * 2.0 * np.pi * t_sec / T_orb

    # Nadir component (rotates once per orbit)
    Q_nad = R3(theta)

    # Spin around velocity axis (local s-hat ≈ y-axis in perifocal frame)
    # Paper Eq. 25: spin matrix around axis aligned with velocity
    Q_spin = np.array([
        [ np.cos(Theta), 0, -np.sin(Theta)],
        [             0, 1,              0],
        [ np.sin(Theta), 0,  np.cos(Theta)],
    ])

    # Eq. 25
    Q_ram = Q_nad @ Q_spin
    # Eq. 26
    Qt = Q.T
    return np.array([Qt @ Q_ram.T @ n for n in N0])


# =============================================================================
# SUN-FIXED ATTITUDES: Sun1, Sun2, Sun3 (Equations 27–28)
# =============================================================================

def attitude_sun_fixed(u_hat, mode='Sun1'):
    """
    Sun-fixed attitude: CubeSat orientation fixed relative to Sun.
    Section 3.3.3, Equations 27–28.

    Parameters
    ----------
    u_hat : ndarray (3,) — unit vector Earth→Sun (from sun_position)
    mode  : str — 'Sun1', 'Sun2', or 'Sun3'

    Returns
    -------
    normals : ndarray (6,3) — face normals in geocentric equatorial frame

    Rotation angles (stated in paper Section 3.3.3):
      Sun1: Theta_s1 = 0°,        Theta_s3 = 0°
      Sun2: Theta_s1 = 0°,        Theta_s3 = 45°
      Sun3: Theta_s1 = 35.2643°,  Theta_s3 = 45°
    """
    Theta_s1_deg = {'Sun1': 0.0,       'Sun2': 0.0,       'Sun3': 35.2643}[mode]
    Theta_s3_deg = {'Sun1': 0.0,       'Sun2': 45.0,      'Sun3': 45.0   }[mode]

    Theta_s1 = np.deg2rad(Theta_s1_deg)
    Theta_s3 = np.deg2rad(Theta_s3_deg)

    u = np.asarray(u_hat, dtype=float)
    u_norm = u / np.linalg.norm(u)

    # Step 1: unit vector perpendicular to Sun in equatorial plane
    # Rotate [u1, u2, 0] by 90° around K-hat
    u_perp = np.array([-u_norm[1], u_norm[0], 0.0])
    perp_norm = np.linalg.norm(u_perp)
    if perp_norm < 1e-9:
        u_perp = np.array([1.0, 0.0, 0.0])
    else:
        u_perp = u_perp / perp_norm

    # Step 2: Eq. 27 — angle between solar vector and K-hat + tilt
    Theta1 = np.arccos(np.clip(u_norm[2], -1, 1)) + Theta_s1

    # Step 3: Rodrigues rotation around u_perp by Theta1
    R_arb = rot_arbitrary(u_perp, Theta1)

    # Step 4: Eq. 28 — angle between u_perp and I-hat (use arctan2 for quadrant)
    Theta2 = np.arctan2(u_perp[1], u_perp[0])

    # Step 5: rotation around K-hat by Theta2
    R_K = R3(-Theta2)

    # Step 6: rotation around J-hat by Theta_s3
    R_J = R2(Theta_s3)

    # Step 7: total rotation (no geocentric transform needed — Section 3.3.3)
    Q_total = R_J @ R_K @ R_arb

    return np.array([Q_total @ n for n in N0])


# =============================================================================
# DISPATCHER
# =============================================================================

def compute_normals(attitude_mode, coe, Q, u_hat, t_sec):
    """
    Return (6,3) array of face normals in geocentric frame for given attitude.

    Parameters
    ----------
    attitude_mode : str   — 'Nadir', 'RAM', 'Sun1', 'Sun2', 'Sun3'
    coe           : dict  — current COE
    Q             : ndarray (3,3) — geocentric→perifocal
    u_hat         : ndarray (3,) — Sun unit vector
    t_sec         : float — elapsed time [s]
    """
    if attitude_mode == 'Nadir':
        return attitude_nadir(coe, Q)
    elif attitude_mode == 'RAM':
        return attitude_ram(coe, Q, t_sec)
    elif attitude_mode in ('Sun1', 'Sun2', 'Sun3'):
        return attitude_sun_fixed(u_hat, mode=attitude_mode)
    else:
        raise ValueError(f"Unknown attitude mode: {attitude_mode}")
