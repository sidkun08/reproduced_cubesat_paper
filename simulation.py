"""
simulation.py
=============
Main orbital propagation loop integrating:
  - Orbit (COE + RK4 + drag + J2)
  - Sun position
  - Shadow detection
  - Attitude normals
  - Irradiance computation
  - Temperature

Reproduces data for Figures 8–12 and Table 2 of:
    Morsch Filho et al. (2020). Energies 13(24), 6691.
"""

import numpy as np
from constants import (MU, RE, DT_SEC, REENTRY_ALT, CD, TLE1, START_DATE,
                       A_1U, A_2U)
from atmosphere import atmospheric_density
from orbital_mechanics import (tle_to_coe, coe_to_position, coe_radius,
                               coe_apogee_perigee, drag_perturbation,
                               j2_perturbation, rk4_step, beta_angle,
                               orbital_velocity_mag)
from sun_shadow import sun_position, check_shadow
from attitude import compute_normals
from irradiance import compute_irradiance, albedo_psi
from thermal import steady_state_temperature


def run_simulation(
    tle=None,
    start_date=None,
    attitude_mode='Sun3',
    dt=DT_SEC,
    max_days=320,
    A_faces=None,
    alpha_over_e=0.5,
    atm_model='USSA76',
    use_drag=True,
    use_j2=True,
    log_every=10,
    A_over_m=0.01654,
):
    """
    Full orbital simulation with irradiance computation.

    Parameters
    ----------
    tle          : dict — TLE parameters (defaults to TLE1 from Table 1)
    start_date   : tuple (year, month, day)
    attitude_mode: str  — 'Nadir','RAM','Sun1','Sun2','Sun3'
    dt           : float — timestep [s]
    max_days     : float — maximum simulation duration [days]
    A_faces      : list  — face areas [m^2], length 6 (default: 1U)
    alpha_over_e : float — absorptivity/emissivity ratio
    atm_model    : str  — 'USSA76' or 'NRLMSISE00'
    use_drag     : bool — include atmospheric drag
    use_j2       : bool — include J2 perturbation
    log_every    : int  — log every N steps (1 = every step)
    A_over_m     : float — frontal area/mass [m^2/kg] for drag

    Returns
    -------
    dict of numpy arrays, one entry per logged quantity.
    """
    if tle is None:
        tle = TLE1
    if start_date is None:
        start_date = START_DATE
    if A_faces is None:
        A_faces = A_1U

    year, month, day = start_date
    coe = tle_to_coe(**tle)
    max_steps = int(max_days * 86400.0 / dt)
    A_total = sum(A_faces)

    logs = {k: [] for k in [
        't_days', 't_s', 'alt', 'alt_apo', 'alt_peri',
        'beta', 'xi', 'e', 'h',
        'Q_sun', 'Q_alb', 'Q_e', 'Q_tot',
        'T_K', 'T_C',
        'Q_N1', 'Q_N2', 'Q_N3', 'Q_N4', 'Q_N5', 'Q_N6',
    ]}

    for step in range(max_steps):
        t = step * dt

        # Sun position (Eq. 21)
        u_hat, r_sun, _, lam, eps = sun_position(year, month, day, UT=t / 3600.0)

        # Satellite position (Eqs. 5–8)
        r_vec, r_mag, Q = coe_to_position(coe)
        alt = r_mag - RE

        if alt < REENTRY_ALT:
            print(f"  [reentry] t={t/86400:.2f} days, alt={alt:.1f} km")
            break

        # Eclipse (Eq. 31)
        xi = check_shadow(r_vec, r_sun)

        # Perturbations (Eqs. 11–18)
        rho = atmospheric_density(alt, model=atm_model)
        prD, psD, pwD = drag_perturbation(coe, rho, CD, A_over_m) if use_drag else (0, 0, 0)
        prJ, psJ, pwJ = j2_perturbation(coe) if use_j2 else (0, 0, 0)
        pr = prD + prJ
        ps = psD + psJ
        pw = pwD + pwJ

        # Integrate (RK4, Eq. 19)
        coe = rk4_step(coe, dt, pr, ps, pw)

        # Log every N steps
        if step % log_every == 0:
            # Attitude normals (Eqs. 22–28)
            normals = compute_normals(attitude_mode, coe, Q, u_hat, t)

            # Irradiance (Eqs. 29–36)
            psi = albedo_psi(r_vec, r_sun)
            irr = compute_irradiance(normals, r_vec, r_sun, xi, A_faces, psi)

            # Temperature (Eq. 38)
            T_K = steady_state_temperature(irr['Q_tot_W'], A_total, alpha_over_e)

            # Beta angle (Eq. 39)
            b = beta_angle(coe, lam, eps)

            # Apogee/perigee
            alt_a, alt_p = coe_apogee_perigee(coe)

            logs['t_days'].append(t / 86400.0)
            logs['t_s'].append(t)
            logs['alt'].append(alt)
            logs['alt_apo'].append(alt_a)
            logs['alt_peri'].append(alt_p)
            logs['beta'].append(b)
            logs['xi'].append(xi)
            logs['e'].append(coe['e'])
            logs['h'].append(coe['h'])
            logs['Q_sun'].append(irr['Q_sun'])
            logs['Q_alb'].append(irr['Q_alb'])
            logs['Q_e'].append(irr['Q_e'])
            logs['Q_tot'].append(irr['Q_tot'])
            logs['T_K'].append(T_K)
            logs['T_C'].append(T_K - 273.15)
            for k in range(6):
                logs[f'Q_N{k+1}'].append(irr['Q_faces_flux'][k])

    return {k: np.array(v) for k, v in logs.items()}


def run_one_orbit(
    tle=None,
    start_date=None,
    t_start_days=0.0,
    attitude_mode='Sun3',
    dt=10.0,
    A_faces=None,
    alpha_over_e=0.5,
    atm_model='NRLMSISE00',
    A_over_m=0.01654,
):
    """
    Simulate exactly one orbital period starting from t_start_days.
    Used for Figures 10–12.

    The simulation first fast-propagates to t_start_days (with dt=60 s),
    then runs one orbital period with fine timestep dt=10 s.

    Returns dict of arrays over one orbit.
    """
    if tle is None:
        tle = TLE1
    if start_date is None:
        start_date = START_DATE
    if A_faces is None:
        A_faces = A_1U

    year, month, day = start_date
    coe = tle_to_coe(**tle)
    A_total = sum(A_faces)

    # --- Fast-propagate to t_start_days ---
    DT_FAST = 60.0
    n_fast   = int(t_start_days * 86400.0 / DT_FAST)
    for step in range(n_fast):
        t = step * DT_FAST
        u_hat, r_sun, _, lam, eps = sun_position(year, month, day, UT=t / 3600.0)
        r_vec, r_mag, Q = coe_to_position(coe)
        alt = r_mag - RE
        if alt < REENTRY_ALT:
            break
        rho = atmospheric_density(alt, model=atm_model)
        prD, psD, pwD = drag_perturbation(coe, rho, CD, A_over_m)
        prJ, psJ, pwJ = j2_perturbation(coe)
        coe = rk4_step(coe, DT_FAST, prD+prJ, psD+psJ, pwD+pwJ)

    # --- Compute orbital period at current state ---
    h, e = coe['h'], coe['e']
    a = h**2 / MU / (1.0 - e**2)
    T_orb = 2.0 * np.pi * np.sqrt(a**3 / MU)   # seconds

    t_offset = t_start_days * 86400.0
    n_orbit  = int(T_orb / dt) + 1

    logs = {k: [] for k in [
        't_s', 'xi',
        'Q_sun', 'Q_alb', 'Q_e', 'Q_tot',
        'T_K', 'T_C',
        'Q_N1', 'Q_N2', 'Q_N3', 'Q_N4', 'Q_N5', 'Q_N6',
    ]}

    coe_orbit = dict(coe)
    for step in range(n_orbit):
        t = t_offset + step * dt
        u_hat, r_sun, _, lam, eps = sun_position(year, month, day, UT=t / 3600.0)
        r_vec, r_mag, Q = coe_to_position(coe_orbit)
        alt = r_mag - RE
        xi = check_shadow(r_vec, r_sun)

        # Propagate
        rho = atmospheric_density(alt, model=atm_model)
        prD, psD, pwD = drag_perturbation(coe_orbit, rho, CD, A_over_m)
        prJ, psJ, pwJ = j2_perturbation(coe_orbit)
        coe_orbit = rk4_step(coe_orbit, dt, prD+prJ, psD+psJ, pwD+pwJ)

        # Irradiance
        normals = compute_normals(attitude_mode, coe_orbit, Q, u_hat, t - t_offset)
        psi = albedo_psi(r_vec, r_sun)
        irr = compute_irradiance(normals, r_vec, r_sun, xi, A_faces, psi)
        T_K = steady_state_temperature(irr['Q_tot_W'], A_total, alpha_over_e)

        logs['t_s'].append(step * dt)
        logs['xi'].append(xi)
        logs['Q_sun'].append(irr['Q_sun'])
        logs['Q_alb'].append(irr['Q_alb'])
        logs['Q_e'].append(irr['Q_e'])
        logs['Q_tot'].append(irr['Q_tot'])
        logs['T_K'].append(T_K)
        logs['T_C'].append(T_K - 273.15)
        for k in range(6):
            logs[f'Q_N{k+1}'].append(irr['Q_faces_flux'][k])

    return {k: np.array(v) for k, v in logs.items()}
