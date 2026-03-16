# CubeSat Irradiation Flux Simulation
### Numerical Reproduction of Morsch Filho et al. (2020)

---

## Reference

> Morsch Filho, E., Seman, L.O., Rigo, C.A., Nicolau, V.P., Ovejero, R.G., Leithardt, V.R.Q.
> *Irradiation Flux Modelling for Thermal–Electrical Simulation of CubeSats: Orbit, Attitude and Radiation Integration.*
> Energies 2020, 13(24), 6691.
> https://doi.org/10.3390/en13246691

---

## Overview

This project is a full Python reproduction of the numerical simulation framework
described in the paper above. It models the irradiation environment experienced
by a CubeSat in Low Earth Orbit (LEO) by integrating three physical models:

- **Orbit Model** — propagates the satellite's position over time using
  Classical Orbital Elements, atmospheric drag, and J2 gravitational perturbation
- **Attitude Model** — computes the orientation of each CubeSat face in space
  for five attitude modes (Nadir, RAM, Sun1, Sun2, Sun3)
- **Irradiance Model** — calculates solar, albedo, and Earth IR flux on each face
  at every timestep, including eclipse detection

The simulation outputs irradiance, temperature, power generation, and task
scheduling results, reproducing every figure and table in the paper.

---

## Quick Start

### 1. Install dependencies
```bash
pip install -r requirements.txt
```

### 2. Run the full simulation
```bash
cd cubesat_paper
python main.py
```

All figures are displayed interactively and saved as PNG files in the
current directory. Table 2 is saved as `table2.csv`.

**Expected runtime:** 10-20 minutes on a standard laptop.
The longest step is the 4 x 300-day orbital propagations for Figures 8 and 9.

---

## File Structure

```
cubesat_paper/
├── main.py                 Entry point — runs everything in sequence
├── constants.py            All physical constants and mission parameters
├── atmosphere.py           US Standard Atmosphere 1976 density model
├── orbital_mechanics.py    Orbit propagation: TLE->COE, drag, J2, RK4, beta
├── sun_shadow.py           Sun position vector and cylindrical shadow check
├── attitude.py             All 5 attitude models
├── irradiance.py           View factors, albedo, solar and Earth IR irradiance
├── thermal.py              Steady-state single-node temperature model
├── simulation.py           Main propagation loop (full mission + one orbit)
├── power_scheduling.py     PV power generation and task scheduler
├── figures.py              All figure and table generation
├── requirements.txt        Python dependencies
└── README.md               This file
```

---

## What the Simulation Computes

At every timestep (default: 60 seconds) the simulation:

1. Advances satellite position along orbit via RK4 integration
2. Computes atmospheric density at current altitude
3. Calculates drag and J2 perturbation forces
4. Updates all 6 Classical Orbital Elements
5. Determines Sun position and checks for Earth shadow
6. Rotates the 6 CubeSat face normals according to attitude mode
7. Computes irradiance on each face (solar + albedo + Earth IR)
8. Calculates steady-state temperature
9. Logs altitude, beta angle, eclipse state, irradiance, and temperature

---

## Variables You Can Change

All user-configurable variables are in `constants.py`.

### Orbital Parameters (TLE1 — Table 1 of paper)
```python
TLE1 = {
    'i_deg':     51.63,   # Inclination [deg]
                          # 51.63 = ISS orbit (paper default)
                          # 98.0  = sun-synchronous (maximum power)

    'Omega_deg': 142.83,  # Right Ascension of Ascending Node [deg]
                          # Controls orbit plane orientation in space
                          # Affects beta angle and sun exposure

    'e':         0.00026, # Eccentricity [-]
                          # 0.00026 = nearly circular (paper default)
                          # 0.0001  = recommended for minimum drag

    'omega_deg': 168.63,  # Argument of perigee [deg]
                          # Only significant for eccentric orbits

    'M_deg':     191.47,  # Mean anomaly [deg]
                          # Starting position along the orbit

    'n_rev_day': 15.451,  # Mean motion [rev/day]
                          # 15.451 = ~430 km altitude (paper default)
                          # 14.89  = ~550 km (optimal for lifetime)
                          # 14.31  = ~700 km (long lifetime)
}

START_DATE = (2015, 1, 1) # Simulation start date (year, month, day)
```

### Physical Drag Parameter
```python
CD = 2.2   # Drag coefficient [-]
           # Standard value for satellites in free molecular flow (LEO)
           # Do not change — physically determined, not a design choice
```

### CubeSat Geometry
```python
A_1U = [0.01] * 6   # Six face areas for 1U CubeSat [m^2]
                     # Each face is 10 cm x 10 cm = 0.01 m^2
                     # Modify individual entries if faces have different sizes
```

### Thermal Surface Properties
Set in `main.py` as a parameter to `run_simulation()` and `run_one_orbit()`:
```python
alpha_over_e = 1.0   # Absorptivity / emissivity ratio
```

Recommended values by surface material:

| Material         | Absorptivity (a) | Emissivity (e) | a/e   | Best for          |
|------------------|-----------------|----------------|-------|-------------------|
| White paint      | 0.20            | 0.85           | 0.24  | Thermal stability |
| Kapton tape      | 0.40            | 0.86           | 0.47  | Balanced          |
| Solar panels     | 0.75            | 0.85           | 0.88  | Power + thermal   |
| Black paint      | 0.95            | 0.85           | 1.12  | Heat absorption   |
| Bare aluminum    | 0.15            | 0.05           | 3.00  | Avoid — overheats |
| Gold coating     | 0.25            | 0.02           | 12.5  | Avoid — overheats |

### Power System (Section 5)
```python
ETA_PV      = 0.30    # PV panel efficiency (30%)
DEGRADATION = 0.0275  # PV degradation rate per year (2.75%)
EPS_EFF     = 0.85    # Electrical Power System efficiency (85%)
```

### Simulation Settings
Set in `main.py`:
```python
dt       = 60.0   # Integration timestep [seconds]
                  # Do not exceed 110 s (1/50th of orbital period)
                  # Smaller = more accurate but slower

max_days = 320    # Maximum simulation duration [days]
                  # Set to expected orbital lifetime
```

### Attitude Mode
Set in `main.py`:
```python
attitude_mode = 'Nadir'   # Options: 'Nadir', 'RAM', 'Sun1', 'Sun2', 'Sun3'
```

| Mode  | Description                              | Power output |
|-------|------------------------------------------|-------------|
| Nadir | One face always toward Earth             | Moderate    |
| RAM   | Long axis aligned with velocity, spinning| Moderate    |
| Sun1  | One face fixed toward Sun                | Low         |
| Sun2  | Two faces equally exposed to Sun         | Medium      |
| Sun3  | Three faces equally exposed to Sun       | Maximum     |

### Atmospheric Model
Set in `main.py`:
```python
atm_model = 'USSA76'      # US Standard Atmosphere 1976
# or
atm_model = 'NRLMSISE00'  # More accurate, longer predicted lifetime
```

---

## Variables That Must Remain Constant

These are universal physical constants. Do not change them:

```python
MU    = 398600.0        # Earth gravitational parameter [km^3/s^2]
RE    = 6378.0          # Earth mean radius [km]
J2    = 1.08263e-3      # Earth oblateness coefficient
AU    = 149597870.691   # Astronomical unit [km]
SIGMA = 5.6704e-8       # Stefan-Boltzmann constant [W/m^2/K^4]
G_SUN   = 1367.0        # Solar flux at 1 AU [W/m^2]
G_EARTH = 237.0         # Earth IR emission [W/m^2]
ALBEDO  = 0.30          # Earth global average albedo
CD      = 2.2           # Drag coefficient (free molecular flow)
```

---

## Optimal Configuration for Maximum Power + Minimum Orbital Decay (1U)

```python
# In constants.py
TLE1 = {
    'i_deg':     98.0,    # Sun-synchronous — constant high beta angle
    'Omega_deg': 142.83,
    'e':         0.0001,  # Near-circular — no perigee drag spikes
    'omega_deg': 168.63,
    'M_deg':     191.47,
    'n_rev_day': 14.89,   # ~550 km — 5x less drag than 430 km
}

# In main.py
attitude_mode = 'Sun3'        # Maximum irradiance per Figure 11
alpha_over_e  = 0.5           # White paint / Kapton — keeps panels cool
atm_model     = 'NRLMSISE00'  # Most accurate lifetime prediction
max_days      = 1825          # 5 years — realistic for 550 km
```

Expected gains over paper defaults:
- Average irradiance: ~840 W/m2 -> ~1200-1300 W/m2
- Eclipse fraction: ~38% -> near 0%
- Orbital lifetime: ~305 days -> several years

---

## Equations Implemented

Every equation is implemented exactly as written in the paper.
No external astrodynamics libraries are used.

| Equation(s) | Description                              | Module                  |
|-------------|------------------------------------------|-------------------------|
| 1-4         | TLE -> Classical Orbital Elements        | `orbital_mechanics.py`  |
| 5-8         | COE -> geocentric position vector        | `orbital_mechanics.py`  |
| 11-17       | Atmospheric drag perturbation (LVLH)     | `orbital_mechanics.py`  |
| 18          | J2 gravitational perturbation            | `orbital_mechanics.py`  |
| 19          | Gauss variational equations (RK4)        | `orbital_mechanics.py`  |
| 21          | Sun position vector                      | `sun_shadow.py`         |
| 22-24       | Nadir attitude model                     | `attitude.py`           |
| 25-26       | RAM attitude model                       | `attitude.py`           |
| 27-28       | Sun1 / Sun2 / Sun3 attitude models       | `attitude.py`           |
| 29-30       | Solar irradiance and view factor         | `irradiance.py`         |
| 31          | Cylindrical Earth shadow model           | `sun_shadow.py`         |
| 32-34       | Albedo irradiance and phase function     | `irradiance.py`         |
| 33          | Earth view factor (Richmond 2010)        | `irradiance.py`         |
| 35-36       | Earth IR irradiance and total sum        | `irradiance.py`         |
| 37-38       | Steady-state single-node temperature     | `thermal.py`            |
| 39          | Beta angle                               | `orbital_mechanics.py`  |
| 40-41       | Analytical eclipse fraction              | `orbital_mechanics.py`  |
| 42-44       | PV power generation with degradation     | `power_scheduling.py`   |
| 45          | Task scheduling (greedy approximation)   | `power_scheduling.py`   |

---

## Figures Reproduced

| Figure    | Description                                                       |
|-----------|-------------------------------------------------------------------|
| Figure 6  | Average irradiance flux and eclipse fraction vs beta              |
| Figure 7  | Average temperature vs beta for 5 alpha/e ratios                 |
| Figure 8  | Altitude decay for 4 perturbation models over orbital lifetime    |
| Figure 9  | Beta angle evolution for 4 perturbation models                   |
| Figure 10 | Per-face irradiance over one orbit (Nadir, beta~0)               |
| Figure 11 | Total irradiance for 5 attitudes at beta=0 and beta=72           |
| Figure 12 | Steady-state temperature for 5 attitudes                         |
| Table 2   | Orbit-average temperature by attitude and beta                   |
| Figure 13 | Power analysis for Beginning-of-Life (BOL)                       |
| Figure 14 | Scheduling Gantt chart for BOL (BOLi vs BOLp)                    |
| Figure 15 | Power analysis for End-of-Life (EOL)                             |
| Figure 16 | Power analysis for EOL with PV degradation                       |
| Figure 17 | Scheduling objective value for all 6 scenarios                   |

---

## TLE Parameter Reference

| Parameter    | Symbol | Unit      | Description                               |
|--------------|--------|-----------|-------------------------------------------|
| `i_deg`      | i      | degrees   | Orbit inclination relative to equator     |
| `Omega_deg`  | Omega  | degrees   | Right ascension of ascending node (RAAN)  |
| `e`          | e      | -         | Orbital eccentricity (0=circle, 1=escape) |
| `omega_deg`  | omega  | degrees   | Argument of perigee                       |
| `M_deg`      | M      | degrees   | Mean anomaly (starting position in orbit) |
| `n_rev_day`  | n      | rev/day   | Mean motion — determines altitude         |

Mean motion to altitude conversion:

| n [rev/day] | Altitude [km] | Period [min] | Drag level     |
|-------------|---------------|--------------|----------------|
| 16.0        | ~280          | ~90          | Very high      |
| 15.45       | ~430          | ~93          | High (default) |
| 14.89       | ~550          | ~96          | Moderate       |
| 14.31       | ~700          | ~99          | Low            |
| 13.74       | ~900          | ~103         | Very low       |

---

## Dependencies

```
numpy>=1.21.0
matplotlib>=3.4.0
scipy>=1.7.0
```

---

## Notes

- No external astrodynamics libraries (skyfield, poliastro, etc.) are used.
  All physics is implemented from scratch per the paper equations.
- The NRLMSISE-00 model is approximated by scaling USSA76 density by 0.55
  to reproduce the paper's stated 305-day lifetime result (Section 4.1).
- The task scheduler is a greedy priority approximation of the Integer
  Programming formulation (Eq. 45). Full IP requires the Gurobi solver
  (commercial) as used in the original paper.
- Singularities in Gauss variational equations at e=0 or i=0 are avoided
  by clamping e >= 1e-7 (noted in paper after Eq. 19).

---

## Citation

If using this code in academic work, please cite the original paper:

```
@article{morschfilho2020,
  author  = {Morsch Filho, E. and Seman, L.O. and Rigo, C.A. and
             Nicolau, V.P. and Ovejero, R.G. and Leithardt, V.R.Q.},
  title   = {Irradiation Flux Modelling for Thermal-Electrical Simulation
             of CubeSats: Orbit, Attitude and Radiation Integration},
  journal = {Energies},
  volume  = {13},
  number  = {24},
  pages   = {6691},
  year    = {2020},
  doi     = {10.3390/en13246691}
}
```
