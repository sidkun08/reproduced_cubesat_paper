"""
figures.py
==========
Generates all figures from:
    Morsch Filho et al. (2020). Energies 13(24), 6691.

Figures reproduced:
  Fig 6  — Irradiance and eclipse fraction vs beta
  Fig 7  — Temperature vs beta for multiple alpha/e ratios
  Fig 8  — Altitude decay for 4 perturbation models
  Fig 9  — Beta angle vs time for 4 perturbation models
  Fig 10 — Per-face irradiance over one orbit (Nadir, beta=0)
  Fig 11 — Total irradiance for 5 attitudes (a: beta=0, b: beta=72)
  Fig 12 — Temperature for 5 attitudes
  Table 2 — Average temperature by attitude and beta
  Fig 13-17 — Power scheduling scenarios
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from constants import TLE1, START_DATE, A_1U, A_2U
from irradiance import orbit_avg_irradiance
from thermal import orbit_avg_temperature, steady_state_temperature
from orbital_mechanics import eclipse_fraction, beta_angle, tle_to_coe
from simulation import run_simulation, run_one_orbit
from power_scheduling import run_all_scenarios, PAYLOADS


ATTITUDES    = ['Nadir', 'RAM', 'Sun1', 'Sun2', 'Sun3']
ATT_COLORS   = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
ATT_LS       = ['-', '-', '-', '-', '-']

ALPHA_E_VALS = [0.5, 1.0, 1.5, 2.0, 2.5]
AE_COLORS    = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

PERT_LABELS  = ['No perturbation', 'Drag:–, J2', 'Drag: USASA76, J2', 'Drag: NRLMSISE-00, J2']
PERT_COLORS  = ['#e377c2', '#2ca02c', '#17becf', '#d62728']

plt.rcParams.update({
    'font.size': 10,
    'axes.grid': True,
    'grid.linestyle': '--',
    'grid.alpha': 0.5,
    'lines.linewidth': 1.5,
    'legend.fontsize': 8,
    'legend.framealpha': 0.9,
})


# =============================================================================
# FIGURE 6 — Irradiance and eclipse fraction vs beta
# =============================================================================

def figure6(save_path=None):
    """
    Reproduce Figure 6: average irradiance flux and eclipse fraction
    as function of beta angle, for altitudes 150 km and 500 km.
    """
    betas = np.linspace(0, 90, 181)
    alts  = [500, 150]

    fig, ax1 = plt.subplots(figsize=(8, 5))
    ax2 = ax1.twinx()

    styles = {500: '-', 150: '--'}
    colors = {
        'sun+alb': '#1f77b4',
        'sun':     '#d62728',
        'e':       '#2ca02c',
        'fE':      '#000000',
    }

    for alt in alts:
        ls = styles[alt]
        Q_sun_arr, Q_alb_arr, Q_e_arr, fE_arr = [], [], [], []
        for b in betas:
            Qs, Qa, Qe, fE = orbit_avg_irradiance(b, alt)
            Q_sun_arr.append(Qs)
            Q_alb_arr.append(Qa)
            Q_e_arr.append(Qe)
            fE_arr.append(fE)
        Q_sun_arr = np.array(Q_sun_arr)
        Q_alb_arr = np.array(Q_alb_arr)

        label_sa = f"Q''_alb+Q''_sun ({alt} km)"
        label_s  = f"Q''_sun ({alt} km)"
        label_e  = f"Q''_e ({alt} km)"
        label_f  = f"fE ({alt} km)"

        ax1.plot(betas, Q_sun_arr + Q_alb_arr, color=colors['sun+alb'], ls=ls, label=label_sa)
        ax1.plot(betas, Q_sun_arr,              color=colors['sun'],     ls=ls, label=label_s)
        ax1.plot(betas, Q_e_arr,                color=colors['e'],       ls=ls, label=label_e)
        ax2.plot(betas, fE_arr,                 color=colors['fE'],      ls=ls, label=label_f)

    ax1.set_xlabel('β [°]')
    ax1.set_ylabel("Irradiance flux [W/m²]")
    ax1.set_xlim(0, 90)
    ax1.set_ylim(0, 1500)
    ax2.set_ylabel("Eclipse fraction fE [–]")
    ax2.set_ylim(0, 0.5)

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2,
               loc='upper left', fontsize=7, ncol=2)

    plt.title('Figure 6 — Irradiance sources and eclipse fraction vs β')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  Saved: {save_path}")
    plt.show()
    plt.close()


# =============================================================================
# FIGURE 7 — Temperature vs beta
# =============================================================================

def figure7(save_path=None):
    """
    Reproduce Figure 7: average steady-state temperature vs beta
    for 5 alpha/e ratios at 150 km and 500 km.
    """
    betas = np.linspace(0, 90, 181)
    alts  = [500, 150]
    # Total external surface area of a 1U CubeSat = 6 faces × 0.01 m^2
    A_total = 6 * 0.01   # m^2

    fig, ax = plt.subplots(figsize=(8, 5))

    for ae_idx, ae in enumerate(ALPHA_E_VALS):
        color = AE_COLORS[ae_idx]
        for alt in alts:
            ls = '-' if alt == 500 else '--'
            T_arr = []
            for b in betas:
                Qs, Qa, Qe, _ = orbit_avg_irradiance(b, alt)
                T = orbit_avg_temperature(Qs, Qa, Qe, A_total, ae, half_area=True)
                T_arr.append(T)
            label = f"α/ε={ae} ({alt} km)" if alt == 500 else None
            ax.plot(betas, T_arr, color=color, ls=ls, label=label)

    ax.set_xlabel('β [°]')
    ax.set_ylabel('Temperature [K]')
    ax.set_xlim(0, 90)
    ax.set_ylim(260, 440)
    ax.legend(loc='upper left', fontsize=7, ncol=2)
    plt.title('Figure 7 — Average temperature vs β for multiple α/ε ratios')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  Saved: {save_path}")
    plt.show()
    plt.close()


# =============================================================================
# FIGURE 8 — Altitude decay for 4 perturbation models
# =============================================================================

def figure8(sim_data_dict, save_path=None):
    """
    Reproduce Figure 8: altitude [km] vs time [days] for 4 perturbation models.

    sim_data_dict: dict mapping PERT_LABELS → simulation result dict
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    for label, color in zip(PERT_LABELS, PERT_COLORS):
        d = sim_data_dict[label]
        t = d['t_days']
        # Show apogee/perigee band for drag cases
        if 'Drag' in label and len(d['alt_apo']) > 0:
            ax.fill_between(t, d['alt_peri'], d['alt_apo'],
                            color=color, alpha=0.6, label=label)
        else:
            ax.plot(t, d['alt'], color=color, lw=1.5, label=label)

    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Altitude (km)')
    ax.set_xlim(0, 310)
    ax.set_ylim(100, 450)
    ax.legend(loc='upper right')
    plt.title('Figure 8 — Altitude evolution for different perturbation models')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  Saved: {save_path}")
    plt.show()
    plt.close()


# =============================================================================
# FIGURE 9 — Beta angle vs time
# =============================================================================

def figure9(sim_data_dict, save_path=None):
    """
    Reproduce Figure 9: beta angle [°] vs time [days] for 4 perturbation models.
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    for label, color in zip(PERT_LABELS, PERT_COLORS):
        d = sim_data_dict[label]
        ax.plot(d['t_days'], d['beta'], color=color, lw=1.2, label=label)

    ax.set_xlabel('Time (days)')
    ax.set_ylabel('β (°)')
    ax.set_xlim(0, 310)
    ax.set_ylim(-100, 100)
    ax.legend(loc='upper right')
    plt.title('Figure 9 — β angle along one year for TLE1')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  Saved: {save_path}")
    plt.show()
    plt.close()


# =============================================================================
# FIGURE 10 — Per-face irradiance over one orbit (Nadir, beta≈0)
# =============================================================================

def figure10(save_path=None):
    """
    Reproduce Figure 10: irradiance on each face vs time for one orbit.
    Nadir attitude, beta=0°, date 10 Jan 2015.
    """
    # t_start=10 days approximates beta≈0 for TLE1 (from Fig 9)
    d = run_one_orbit(attitude_mode='Nadir', t_start_days=10.0,
                      A_faces=A_1U, alpha_over_e=1.0)

    face_colors = ['#1f77b4','#ff7f0e','#2ca02c','#d62728','#9467bd','#8c564b']
    fig, ax = plt.subplots(figsize=(8, 5))

    for k in range(6):
        ax.plot(d['t_s'], d[f'Q_N{k+1}'], color=face_colors[k],
                lw=1.3, label=f'N{k+1}')

    ax.set_xlabel('Time [s]')
    ax.set_ylabel("Q''_sun + Q''_alb + Q''_e  [W/m²]")
    ax.set_xlim(0, d['t_s'][-1])
    ax.set_ylim(0, 1500)
    ax.legend(loc='upper left', ncol=2)
    plt.title('Figure 10 — Per-face irradiance (Nadir, β≈0°, 10 Jan 2015)')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  Saved: {save_path}")
    plt.show()
    plt.close()


# =============================================================================
# FIGURE 11 — Total irradiance for 5 attitudes (2 subplots)
# =============================================================================

def figure11(save_path=None):
    """
    Reproduce Figure 11: total irradiance for 5 attitudes.
    (a) beta=0° at 10 days,  (b) beta=72° at 143 days.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

    configs = [
        (10.0,  'β = 0°, at 10 days'),
        (143.0, 'β = 72°, at 143 days'),
    ]

    for ax, (t_days, subtitle) in zip(axes, configs):
        for att, color in zip(ATTITUDES, ATT_COLORS):
            d = run_one_orbit(attitude_mode=att, t_start_days=t_days,
                              A_faces=A_1U, alpha_over_e=1.0)
            ax.plot(d['t_s'], d['Q_tot'], color=color, lw=1.3, label=att)

        ax.set_xlabel('Time [s]')
        ax.set_xlim(0, 5500)
        ax.set_ylim(0, 4000)
        ax.set_title(f'({subtitle})')
        ax.legend(loc='upper left')

    axes[0].set_ylabel("Q''_tot [W/m²]")
    fig.suptitle('Figure 11 — Total irradiance flux for each attitude')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  Saved: {save_path}")
    plt.show()
    plt.close()


# =============================================================================
# FIGURE 12 — Temperature for 5 attitudes (2 subplots)
# =============================================================================

def figure12(save_path=None):
    """
    Reproduce Figure 12: steady-state temperature for 5 attitudes.
    α/ε = 1.0 (stated in Section 4.1).
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

    configs = [
        (10.0,  'β = 0°, at 10 days'),
        (143.0, 'β = 72°, at 143 days'),
    ]

    for ax, (t_days, subtitle) in zip(axes, configs):
        for att, color in zip(ATTITUDES, ATT_COLORS):
            d = run_one_orbit(attitude_mode=att, t_start_days=t_days,
                              A_faces=A_1U, alpha_over_e=1.0)
            ax.plot(d['t_s'], d['T_K'], color=color, lw=1.3, label=att)

        ax.set_xlabel('Time [s]')
        ax.set_xlim(0, 5500)
        ax.set_ylim(180, 340)
        ax.set_title(f'({subtitle})')
        ax.legend(loc='upper left')

    axes[0].set_ylabel('Temperature [K]')
    fig.suptitle('Figure 12 — Steady-state temperature for each attitude (α/ε=1.0)')
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"  Saved: {save_path}")
    plt.show()
    plt.close()


# =============================================================================
# TABLE 2 — Average temperature by attitude and beta
# =============================================================================

def table2(save_path=None):
    """
    Reproduce Table 2: orbit-average temperature for each attitude and beta.
    """
    print("\n" + "=" * 65)
    print("  TABLE 2 — Average Temperature [K] for each Attitude and β")
    print("=" * 65)
    print(f"  {'β':>5} | {'Nadir':>8} {'RAM':>8} {'Sun1':>8} {'Sun2':>8} {'Sun3':>8}")
    print("  " + "-" * 55)

    results = {}
    for t_days, beta_label in [(10.0, 0), (143.0, 72)]:
        row = {}
        for att in ATTITUDES:
            d = run_one_orbit(attitude_mode=att, t_start_days=t_days,
                              A_faces=A_1U, alpha_over_e=1.0)
            avg_T = np.mean(d['T_K'])
            row[att] = avg_T
        results[beta_label] = row
        vals = [f"{row[a]:>8.0f}" for a in ATTITUDES]
        print(f"  {beta_label:>5} | {''.join(vals)}")

    print("=" * 65)
    print("  Paper values:")
    print("       0 |     256      258      248      259      265")
    print("      72 |     289      293      274      292      304")
    print("=" * 65)

    if save_path:
        with open(save_path, 'w') as f:
            f.write("beta,Nadir,RAM,Sun1,Sun2,Sun3\n")
            for beta_label, row in results.items():
                vals = ','.join(f"{row[a]:.1f}" for a in ATTITUDES)
                f.write(f"{beta_label},{vals}\n")
        print(f"  Saved: {save_path}")

    return results


# =============================================================================
# FIGURES 13-17 — Power scheduling
# =============================================================================

def figures13_to_17(save_prefix=None):
    """
    Reproduce Figures 13–17: power profiles and scheduling for 6 scenarios.
    """
    print("\n  Running all 6 scheduling scenarios (Section 5)...")
    scenarios = run_all_scenarios()

    scenario_pairs = [
        ('BOLi', 'BOLp',   'Figure 13 — Power analysis for BOL'),
        ('EOLi', 'EOLp',   'Figure 15 — Power analysis for EOL'),
        ('EOLi,d','EOLp,d','Figure 16 — Power analysis for EOL with degradation'),
    ]

    for pair_idx, (sc_ideal, sc_prop, title) in enumerate(scenario_pairs):
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))

        for ax, sc_name in zip(axes, [sc_ideal, sc_prop]):
            sc = scenarios[sc_name]
            t_min = np.arange(len(sc['power']))
            ax.step(t_min, sc['power'], where='post', color='#1f77b4',
                    lw=1.2, label='Available power', linestyle=':')
            ax.step(t_min, sc['usage'], where='post', color='#d62728',
                    lw=1.5, label="Payload's total usage")
            ax.set_xlabel('Time [min]')
            ax.set_ylabel('Power [W]')
            ax.set_title(f'({sc_name})')
            ax.set_xlim(0, len(sc['power']))
            ax.set_ylim(0, 10)
            ax.legend(loc='upper right', fontsize=8)

        fig.suptitle(title)
        plt.tight_layout()
        if save_prefix:
            path = f"{save_prefix}_fig{13+pair_idx}.png"
            plt.savefig(path, dpi=150)
            print(f"  Saved: {path}")
        plt.show()
        plt.close()

    # Figure 14: scheduling Gantt for BOL
    fig14, axes14 = plt.subplots(1, 2, figsize=(13, 6))
    for ax, sc_name in zip(axes14, ['BOLi', 'BOLp']):
        sc = scenarios[sc_name]
        t_len = len(sc['power'])
        payload_ids = [p[0] for p in PAYLOADS]
        for y_idx, pid in enumerate(payload_ids):
            for (start, dur) in sc['schedule'].get(pid, []):
                ax.barh(y_idx, dur, left=start, height=0.6,
                        color='#1f77b4', edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(payload_ids)))
        ax.set_yticklabels(payload_ids)
        ax.set_xlabel('Time [min]')
        ax.set_xlim(0, t_len)
        ax.set_title(f'Optimal Scheduling ({sc_name})')
        ax.set_ylabel('On/Off')
    fig14.suptitle('Figure 14 — Scheduling results for BOL')
    plt.tight_layout()
    if save_prefix:
        path = f"{save_prefix}_fig14.png"
        plt.savefig(path, dpi=150)
        print(f"  Saved: {path}")
    plt.show()
    plt.close()

    # Figure 17: objective value bar chart
    sc_names  = list(scenarios.keys())
    obj_vals  = [scenarios[s]['obj'] for s in sc_names]

    fig17, ax17 = plt.subplots(figsize=(9, 4))
    bars = ax17.bar(sc_names, obj_vals, color='#2ca02c', edgecolor='black', width=0.5)
    for bar, val in zip(bars, obj_vals):
        ax17.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                  f'{val:.0f}', ha='center', va='bottom', fontsize=8)
    ax17.set_ylabel('Objective value')
    ax17.set_ylim(500, max(obj_vals) * 1.15)
    ax17.set_title('Figure 17 — Scheduling objective value for each scenario')
    plt.tight_layout()
    if save_prefix:
        path = f"{save_prefix}_fig17.png"
        plt.savefig(path, dpi=150)
        print(f"  Saved: {path}")
    plt.show()
    plt.close()
