"""
constants.py
============
Physical and mission constants for reproduction of:
    Morsch Filho et al. (2020). Energies 13(24), 6691.
    https://doi.org/10.3390/en13246691

All values taken verbatim from paper text/equations.
"""

# --- Orbital mechanics (Section 3, paper) ---
MU      = 398600.0       # km^3/s^2  — Earth gravitational parameter
RE      = 6378.0         # km        — Earth mean radius
J2      = 1.08263e-3     # [-]       — Second zonal harmonic (Eq. 18)

# --- Radiation sources (Section 2, paper) ---
G_SUN   = 1367.0         # W/m^2    — Solar flux at 1 AU
G_EARTH = 237.0          # W/m^2    — Earth IR flux
ALBEDO  = 0.30           # [-]      — Global annual albedo coefficient

# --- Astronomical (Eq. 21h) ---
AU      = 149597870.691  # km       — Astronomical unit

# --- Thermal (Eq. 37-38) ---
SIGMA   = 5.6704e-8      # W/m^2/K^4 — Stefan-Boltzmann constant

# --- Drag (Section 3.1.1) ---
CD      = 2.2            # [-]      — Drag coefficient (constant, per Haneveer 2017)

# --- TLE1: PropCube-2, Table 1, paper ---
TLE1 = {
    'i_deg':     98.0, # original: 51.63
    'Omega_deg': 150.0, # Original: 142.83
    'e':         0.00026,
    'omega_deg': 168.63,
    'M_deg':     191.47,
    'n_rev_day': 15.451,
}
START_DATE = (2015, 1, 1)  # year, month, day — Section 4

# --- CubeSat geometry ---
# 1U: all faces 10x10 cm (Section 3, "without deployable parts")
A_1U = [0.01] * 6   # m^2, six faces

# 2U: faces from Section 5
A_2U = [0.012072, 0.012072, 0.006036, 0.006036, 0.012072, 0.012072]  # m^2

# --- Power/PV (Section 5) ---
ETA_PV       = 0.30    # PV efficiency
DEGRADATION  = 0.0275  # per year (2.75%)
EPS_EFF      = 0.85    # Electrical Power System efficiency

# --- Attitude spin rate (Section 4, RAM case) ---
RAM_SPIN_REV_PER_ORBIT = 4  # rotations per orbit

# --- Simulation (Section 4) ---
DT_SEC      = 60.0   # integration timestep [s]
REENTRY_ALT = 100.0  # km — reentry threshold
