"""
main.py
=======
Entry point for full reproduction of:
    Morsch Filho et al. (2020). Energies 13(24), 6691.
    https://doi.org/10.3390/en13246691

Run:
    python main.py

Generates all figures and Table 2. Figures are displayed interactively
and saved as PNG files in the current directory.

Execution time: approximately 10–20 minutes on a standard laptop.
"""

import os
import sys
import time
import numpy as np

# Ensure local modules are found
sys.path.insert(0, os.path.dirname(__file__))

from constants import TLE1, START_DATE
from simulation import run_simulation
from figures import (figure6, figure7, figure8, figure9,
                     figure10, figure11, figure12, table2,
                     figures13_to_17,
                     PERT_LABELS)


def validate_initial_conditions():
    """
    Validate key computed quantities against paper-stated values.
    Prints pass/fail for each checkpoint.
    """
    from orbital_mechanics import tle_to_coe, coe_to_position, coe_radius
    from atmosphere import atmospheric_density
    from sun_shadow import sun_position

    print("\n" + "=" * 60)
    print("  VALIDATION CHECKPOINTS")
    print("=" * 60)

    coe = tle_to_coe(**TLE1)

    # 1. Initial altitude
    r_vec, r_mag, Q = coe_to_position(coe)
    alt = r_mag - 6378.0
    status = "PASS" if 420 < alt < 440 else "FAIL"
    print(f"  [{status}] Initial altitude: {alt:.1f} km  (expected ~430 km)")

    # 2. Orbital period
    n = TLE1['n_rev_day'] * 2 * np.pi / 86400.0
    T = 2 * np.pi / n
    status = "PASS" if 5500 < T < 5700 else "FAIL"
    print(f"  [{status}] Orbital period:   {T:.0f} s  (expected ~5572 s)")

    # 3. Eccentricity handled
    status = "PASS" if coe['e'] > 1e-7 else "FAIL"
    print(f"  [{status}] Eccentricity:     {coe['e']:.5f}  (expected 0.00026)")

    # 4. Atmosphere at 430 km
    rho = atmospheric_density(alt, model='USSA76')
    status = "PASS" if 1e-12 < rho < 1e-11 else "FAIL"
    print(f"  [{status}] USSA76 @ {alt:.0f} km: {rho:.3e} kg/m³  (expected ~3.7e-12)")

    # 5. Sun position Jan 1 2015
    u_hat, _, _, lam, eps = sun_position(2015, 1, 1, UT=0.0)
    lam_deg = np.rad2deg(lam) % 360
    status = "PASS" if 270 < lam_deg < 290 else "FAIL"
    print(f"  [{status}] Sun ecliptic lon: {lam_deg:.1f}°  (expected ~280° for Jan 1)")

    # 6. Eclipse fraction at beta=0, 430 km
    from orbital_mechanics import eclipse_fraction
    fE = eclipse_fraction(0.0, alt)
    status = "PASS" if 0.30 < fE < 0.40 else "FAIL"
    print(f"  [{status}] Eclipse frac β=0: {fE:.3f}  (expected ~0.35)")

    print("=" * 60 + "\n")


def run_perturbation_simulations():
    """
    Run all 4 perturbation model simulations needed for Figures 8 and 9.
    Returns dict: label → simulation result.

    Paper Section 4.1 states:
      - USSA76 + J2: reentry ~196 days
      - NRLMSISE-00 + J2: reentry ~305 days
    """
    print("  Running 4 perturbation simulations (Figures 8 & 9)...")
    print("  [This is the most time-consuming step — ~10-15 min total]\n")

    configs = [
        ('No perturbation',         dict(use_drag=False, use_j2=False, atm_model='USSA76')),
        ('Drag:–, J2',              dict(use_drag=False, use_j2=True,  atm_model='USSA76')),
        ('Drag: USASA76, J2',       dict(use_drag=True,  use_j2=True,  atm_model='USSA76')),
        ('Drag: NRLMSISE-00, J2',   dict(use_drag=True,  use_j2=True,  atm_model='NRLMSISE00')),
    ]

    results = {}
    for label, kwargs in configs:
        t0 = time.time()
        print(f"  Simulating: {label} ...", end='', flush=True)
        d = run_simulation(
            tle=TLE1,
            start_date=START_DATE,
            attitude_mode='Sun3',
            dt=60.0,
            max_days=320,
            log_every=10,
            **kwargs,
        )
        elapsed = time.time() - t0
        final_alt = d['alt'][-1] if len(d['alt']) > 0 else 0
        lifetime  = d['t_days'][-1] if len(d['t_days']) > 0 else 0
        print(f"  done ({elapsed:.0f}s). Lifetime: {lifetime:.0f} days, final alt: {final_alt:.0f} km")
        results[label] = d

    return results


def main():
    print("=" * 60)
    print("  CubeSat Irradiation Simulation — Paper Reproduction")
    print("  Morsch Filho et al. (2020), Energies 13(24), 6691")
    print("=" * 60)

    # Step 0: Validate
    validate_initial_conditions()

    # Step 1: Figures 6 & 7 (analytical, fast)
    print("  Generating Figure 6 (irradiance vs beta)...")
    figure6(save_path='figure6.png')

    print("  Generating Figure 7 (temperature vs beta)...")
    figure7(save_path='figure7.png')

    # Step 2: Run perturbation simulations (slow — ~10-15 min)
    sim_data = run_perturbation_simulations()

    # Step 3: Figures 8 & 9
    print("  Generating Figure 8 (altitude decay)...")
    figure8(sim_data, save_path='figure8.png')

    print("  Generating Figure 9 (beta angle vs time)...")
    figure9(sim_data, save_path='figure9.png')

    # Step 4: One-orbit figures (Figures 10-12)
    print("  Generating Figure 10 (per-face irradiance, Nadir)...")
    figure10(save_path='figure10.png')

    print("  Generating Figure 11 (total irradiance, 5 attitudes)...")
    figure11(save_path='figure11.png')

    print("  Generating Figure 12 (temperature, 5 attitudes)...")
    figure12(save_path='figure12.png')

    # Step 5: Table 2
    print("  Computing Table 2 (average temperature)...")
    table2(save_path='table2.csv')

    # Step 6: Scheduling figures (13-17)
    print("  Generating Figures 13-17 (power scheduling)...")
    figures13_to_17(save_prefix='figure')

    print("\n" + "=" * 60)
    print("  All figures generated successfully.")
    print("  Files saved in current directory.")
    print("=" * 60)


if __name__ == '__main__':
    main()
