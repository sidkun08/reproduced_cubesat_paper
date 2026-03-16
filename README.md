# CubeSat Irradiation Flux Simulation
## Reproduction of Morsch Filho et al. (2020), *Energies* 13(24), 6691

### Reference
> Morsch Filho, E., Seman, L.O., Rigo, C.A., Nicolau, V.P., Ovejero, R.G., Leithardt, V.R.Q.  
> *Irradiation Flux Modelling for Thermal–Electrical Simulation of CubeSats: Orbit, Attitude and Radiation Integration.*  
> Energies 2020, 13(24), 6691. https://doi.org/10.3390/en13246691

---

### File Structure

| File | Contents |
|---|---|
| `constants.py` | All physical constants and TLE1 data from paper |
| `atmosphere.py` | US Standard Atmosphere 1976 density model |
| `orbital_mechanics.py` | TLE→COE, COE→position, drag, J2, Gauss eqs, RK4, beta angle |
| `sun_shadow.py` | Sun position (Eq. 21), cylindrical shadow (Eq. 31) |
| `attitude.py` | Nadir, RAM, Sun1, Sun2, Sun3 attitude models (Eqs. 22–28) |
| `irradiance.py` | View factors, albedo psi, irradiance (Eqs. 29–36) |
| `thermal.py` | Steady-state temperature (Eqs. 37–38) |
| `simulation.py` | Main propagation loop |
| `power_scheduling.py` | PV power (Eqs. 42–44), scheduler, 6 scenarios (Table 4) |
| `figures.py` | All figure generation (Figs 6–17, Table 2) |
| `main.py` | Entry point — runs all figures in sequence |

---

### Installation

```bash
pip install -r requirements.txt
```

### Usage

```bash
cd cubesat_paper
python main.py
```

All figures are displayed and saved as PNG files in the current directory.

**Expected runtime:** 10–20 minutes (dominated by the 4 × 300-day orbital simulations for Figures 8 & 9).

---

### Equations Implemented

| Equation(s) | Description | Module |
|---|---|---|
| 1–4 | TLE → Classical Orbital Elements | `orbital_mechanics.py` |
| 5–8 | COE → geocentric position | `orbital_mechanics.py` |
| 11–17 | Atmospheric drag perturbation (LVLH) | `orbital_mechanics.py` |
| 18 | J2 gravitational perturbation | `orbital_mechanics.py` |
| 19 | Gauss variational equations (RK4) | `orbital_mechanics.py` |
| 21 | Sun position vector | `sun_shadow.py` |
| 22–24 | Nadir attitude | `attitude.py` |
| 25–26 | RAM attitude | `attitude.py` |
| 27–28 | Sun1/Sun2/Sun3 attitudes | `attitude.py` |
| 29, 30 | Solar irradiance + view factor | `irradiance.py` |
| 31 | Cylindrical shadow model | `sun_shadow.py` |
| 32–34 | Albedo irradiance + phase function | `irradiance.py` |
| 33 | Earth view factor (Richmond 2010) | `irradiance.py` |
| 35–36 | Earth IR + total irradiance | `irradiance.py` |
| 37–38 | Steady-state temperature | `thermal.py` |
| 39 | Beta angle | `orbital_mechanics.py` |
| 40–41 | Eclipse fraction | `orbital_mechanics.py` |
| 42–44 | PV power + degradation | `power_scheduling.py` |

---

### Figures Reproduced

- **Figure 6:** Average irradiance flux and eclipse fraction vs β for 150 km and 500 km
- **Figure 7:** Average temperature vs β for 5 α/ε ratios at two altitudes
- **Figure 8:** Altitude decay for 4 perturbation models over orbital lifetime
- **Figure 9:** β angle evolution for 4 perturbation models
- **Figure 10:** Per-face irradiance over one orbit (Nadir, β≈0°)
- **Figure 11:** Total irradiance for 5 attitude modes (β=0° and β=72°)
- **Figure 12:** Steady-state temperature for 5 attitude modes
- **Table 2:** Orbit-average temperature by attitude and β
- **Figures 13–17:** Power scheduling scenarios (BOL/EOL, ideal/proposed)
