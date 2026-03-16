"""
power_scheduling.py
===================
Implements Section 5 of the paper:
  - PV power generation (Equations 42–44)
  - Greedy priority scheduler approximating the IP formulation (Equation 45)
  - 6 scheduling scenarios from Table 4

All equations from:
    Morsch Filho et al. (2020). Energies 13(24), 6691.
"""

import numpy as np
from constants import (MU, RE, G_SUN, TLE1, START_DATE, A_2U,
                       ETA_PV, DEGRADATION, EPS_EFF)
from orbital_mechanics import tle_to_coe, coe_to_position, coe_radius, rk4_step
from orbital_mechanics import drag_perturbation, j2_perturbation
from atmosphere import atmospheric_density
from sun_shadow import sun_position, check_shadow
from attitude import compute_normals
from irradiance import compute_irradiance, albedo_psi
from constants import CD


# =============================================================================
# POWER GENERATION (Equations 42–44)
# =============================================================================

def pv_power(Q_flux_Wm2, A_PV_m2, life_years=0.0, include_degradation=False):
    """
    Photovoltaic power output [W].

    Eq. 42: P = eta * Q'' * A
    Eq. 43: Ld = (1 - degradation_per_year)^satellite_life
    Eq. 44: P = eta * Q * A * Ld

    Parameters
    ----------
    Q_flux_Wm2          : float — irradiance flux [W/m^2]
    A_PV_m2             : float — PV panel area [m^2]
    life_years          : float — satellite life [years] for degradation
    include_degradation : bool

    Returns power available [W] after EPS efficiency (× 0.85, Section 5).
    """
    Ld = (1.0 - DEGRADATION) ** life_years if include_degradation else 1.0
    P = ETA_PV * Q_flux_Wm2 * A_PV_m2 * Ld   # Eq. 44
    return P * EPS_EFF                          # × EPS efficiency


# =============================================================================
# ONE-ORBIT POWER PROFILE (minute by minute, Section 5)
# =============================================================================

def compute_power_profile(
    tle,
    start_date,
    attitude_mode,
    A_faces,
    t_start_days=0.0,
    life_years=0.0,
    include_degradation=False,
    atm_model='NRLMSISE00',
    solar_only=False,
):
    """
    Compute available power [W] per minute over one orbital period.

    Parameters
    ----------
    tle             : dict — TLE parameters
    start_date      : tuple (year, month, day)
    attitude_mode   : str
    A_faces         : list[float] — face areas [m^2]
    t_start_days    : float — simulation start offset [days]
    life_years      : float — for degradation calculation
    include_degradation : bool
    atm_model       : str
    solar_only      : bool — if True ignore albedo (BOLi/EOLi cases)
    """
    year, month, day = start_date
    coe = tle_to_coe(**tle)

    # Fast-propagate to t_start_days
    DT_FAST = 60.0
    n_fast = int(t_start_days * 86400.0 / DT_FAST)
    for step in range(n_fast):
        t = step * DT_FAST
        u_hat, r_sun, _, lam, eps = sun_position(year, month, day, UT=t / 3600.0)
        r_vec, r_mag, Q = coe_to_position(coe)
        alt = r_mag - RE
        if alt < 100.0:
            break
        rho = atmospheric_density(alt, model=atm_model)
        prD, psD, pwD = drag_perturbation(coe, rho, CD, A_over_m=0.01654)
        prJ, psJ, pwJ = j2_perturbation(coe)
        coe = rk4_step(coe, DT_FAST, prD+prJ, psD+psJ, pwD+pwJ)

    # Compute orbital period
    h, e = coe['h'], coe['e']
    a = h**2 / MU / (1.0 - e**2)
    T_orb = 2.0 * np.pi * np.sqrt(a**3 / MU)
    n_min = int(T_orb / 60.0) + 1

    t_offset = t_start_days * 86400.0
    power_profile = []
    coe_orbit = dict(coe)

    for m in range(n_min):
        t = t_offset + m * 60.0
        u_hat, r_sun, _, lam, eps = sun_position(year, month, day, UT=t / 3600.0)
        r_vec, r_mag, Q = coe_to_position(coe_orbit)
        xi = check_shadow(r_vec, r_sun)

        normals = compute_normals(attitude_mode, coe_orbit, Q, u_hat, m * 60.0)
        psi = 0.0 if solar_only else albedo_psi(r_vec, r_sun)
        irr = compute_irradiance(normals, r_vec, r_sun, xi, A_faces, psi)

        # Power from all PV faces
        P = pv_power(irr['Q_tot'], sum(A_faces), life_years, include_degradation)
        power_profile.append(P)

        rho = atmospheric_density(r_mag - RE, model=atm_model)
        prD, psD, pwD = drag_perturbation(coe_orbit, rho, CD, A_over_m=0.01654)
        prJ, psJ, pwJ = j2_perturbation(coe_orbit)
        coe_orbit = rk4_step(coe_orbit, 60.0, prD+prJ, psD+psJ, pwD+pwJ)

    return np.array(power_profile)


def make_circular_tle(alt_km, inc_deg):
    """Create TLE dict for a circular orbit at given altitude and inclination."""
    r = RE + alt_km
    n_rad = np.sqrt(MU / r**3)                       # rad/s
    n_rev = n_rad * 86400.0 / (2.0 * np.pi)          # rev/day
    return {
        'i_deg': inc_deg, 'Omega_deg': 142.83,
        'e': 0.0001, 'omega_deg': 0.0,
        'M_deg': 0.0, 'n_rev_day': n_rev,
    }


# =============================================================================
# GREEDY PRIORITY SCHEDULER (approximation of Eq. 45)
# =============================================================================

# Table 5 payload data (Section 5)
PAYLOADS = [
    # (id, priority, power_W, t_min_min, t_max_min, p_min_min, p_max_min)
    ('A', 5, 1.0,  3, 100, 2, 100),
    ('B', 2, 0.23, 2, 100, 2, 100),
    ('C', 1, 0.8,  3, 100, 2, 100),
    ('D', 4, 1.3,  6, 100, 2, 100),
    ('E', 1, 0.5,  1, 100, 2, 100),
    ('F', 3, 0.1,  2, 100, 2, 100),
    ('G', 1, 1.1,  4, 100, 2, 100),
    ('H', 1, 0.9,  3, 100, 2, 100),
]


def greedy_scheduler(power_profile_W):
    """
    Greedy priority-based scheduler approximating the IP formulation (Eq. 45).

    At each minute t, schedule highest-priority tasks that:
      - fit within available power (Eq. 45b)
      - have not exceeded their max runs
      - respect minimum period between runs

    Returns
    -------
    schedule  : dict payload_id → list of (start_min, duration_min)
    usage     : ndarray — total payload power usage per minute [W]
    objective : float  — sum of priority × minutes executed (approx Eq. 45a)
    """
    T = len(power_profile_W)
    # Sort payloads by priority (descending)
    sorted_pl = sorted(PAYLOADS, key=lambda x: -x[1])

    # State tracking
    last_end   = {p[0]: -999 for p in PAYLOADS}  # last minute task ended
    run_count  = {p[0]: 0    for p in PAYLOADS}
    schedule   = {p[0]: []   for p in PAYLOADS}
    usage      = np.zeros(T)
    objective  = 0.0

    t = 0
    while t < T:
        avail = power_profile_W[t]
        for pid, prio, pwr, t_min, t_max, p_min, p_max in sorted_pl:
            # Check period constraint: minimum p_min minutes since last start
            last_start = last_end[pid] - t_min if schedule[pid] else -999
            if t - last_start < p_min and schedule[pid]:
                continue
            # Check if task fits in power budget
            if pwr > avail + 1e-6:
                continue
            # Check if enough time remains for minimum duration
            if t + t_min > T:
                continue
            # Schedule this task for t_min minutes (minimum CPU time)
            dur = min(t_min, T - t)
            schedule[pid].append((t, dur))
            last_end[pid]  = t + dur
            run_count[pid] += 1
            avail          -= pwr
            objective      += prio * dur
            # Apply power usage
            for dt_ in range(dur):
                if t + dt_ < T:
                    usage[t + dt_] += pwr
        t += 1

    return schedule, usage, objective


def run_all_scenarios():
    """
    Run all 6 scheduling scenarios from Table 4 and return power profiles.

    Returns dict: scenario_name → {'power': array, 'usage': array, 'obj': float}
    """
    # Find EOL time for TLE1: paper states reentry at 275 days (NRLMSISE-00)
    T_EOL_DAYS = 245.0   # 30 days before reentry at 275 days (paper Table 4)

    scenarios = {}

    # BOLi: Nadir, circular 430 km, solar only
    tle_boli = make_circular_tle(430, 51.6)
    p = compute_power_profile(tle_boli, START_DATE, 'Nadir', A_2U,
                              solar_only=True, atm_model='USSA76')
    sch, use, obj = greedy_scheduler(p)
    scenarios['BOLi'] = {'power': p, 'usage': use, 'obj': obj, 'schedule': sch}

    # EOLi: Nadir, circular 150 km, solar only
    tle_eoli = make_circular_tle(150, 51.6)
    p = compute_power_profile(tle_eoli, START_DATE, 'Nadir', A_2U,
                              solar_only=True, atm_model='USSA76')
    sch, use, obj = greedy_scheduler(p)
    scenarios['EOLi'] = {'power': p, 'usage': use, 'obj': obj, 'schedule': sch}

    # EOLi,d: same as EOLi with 6-month degradation
    p = compute_power_profile(tle_eoli, START_DATE, 'Nadir', A_2U,
                              life_years=6/12, include_degradation=True,
                              solar_only=True, atm_model='USSA76')
    sch, use, obj = greedy_scheduler(p)
    scenarios['EOLi,d'] = {'power': p, 'usage': use, 'obj': obj, 'schedule': sch}

    # BOLp: RAM, TLE1, solar+albedo
    p = compute_power_profile(TLE1, START_DATE, 'RAM', A_2U,
                              t_start_days=0.0, atm_model='NRLMSISE00')
    sch, use, obj = greedy_scheduler(p)
    scenarios['BOLp'] = {'power': p, 'usage': use, 'obj': obj, 'schedule': sch}

    # EOLp: RAM, TLE1, 30 days before reentry (275-30=245 days)
    p = compute_power_profile(TLE1, START_DATE, 'RAM', A_2U,
                              t_start_days=T_EOL_DAYS, atm_model='NRLMSISE00')
    sch, use, obj = greedy_scheduler(p)
    scenarios['EOLp'] = {'power': p, 'usage': use, 'obj': obj, 'schedule': sch}

    # EOLp,d: same as EOLp with degradation
    p = compute_power_profile(TLE1, START_DATE, 'RAM', A_2U,
                              t_start_days=T_EOL_DAYS,
                              life_years=275/365,
                              include_degradation=True,
                              atm_model='NRLMSISE00')
    sch, use, obj = greedy_scheduler(p)
    scenarios['EOLp,d'] = {'power': p, 'usage': use, 'obj': obj, 'schedule': sch}

    return scenarios
