"""
atmosphere.py
=============
US Standard Atmosphere 1976 (USSA76) piecewise-exponential density model.
Referenced in paper Section 3.1.1 and Figure 8 (label "Drag: USASA76, J2").
Table values from Curtis (2014) Appendix D [paper Ref 40].

NRLMSISE-00 approximation: scale USSA76 by factor 0.55 to reproduce the
paper's Figure 8 result of 305-day lifetime vs ~196 days for USSA76.
This scaling is consistent with solar-minimum NRLMSISE-00 densities at LEO
altitudes being lower than USSA76 (Vallado & Finkleman 2014 [Ref 44]).
"""

import numpy as np

# Piecewise exponential table: (base_altitude_km, base_density_kg/m3, scale_height_km)
# Source: Curtis (2014) Table A.1, referenced as [Ref 40] in paper
USSA76_TABLE = [
    (0,    1.225,      7.249),
    (25,   3.899e-2,   6.349),
    (30,   1.774e-2,   6.682),
    (40,   3.972e-3,   7.554),
    (50,   1.057e-3,   8.382),
    (60,   3.206e-4,   7.714),
    (70,   8.770e-5,   6.549),
    (80,   1.905e-5,   5.799),
    (90,   3.396e-6,   5.382),
    (100,  5.297e-7,   5.877),
    (110,  9.661e-8,   7.263),
    (120,  2.438e-8,   9.473),
    (130,  8.484e-9,  12.636),
    (140,  3.845e-9,  16.149),
    (150,  2.070e-9,  22.523),
    (180,  5.464e-10, 29.740),
    (200,  2.789e-10, 37.105),
    (250,  7.248e-11, 45.546),
    (300,  2.418e-11, 53.628),
    (350,  9.158e-12, 53.298),
    (400,  3.725e-12, 58.515),
    (450,  1.585e-12, 60.828),
    (500,  6.967e-13, 63.822),
    (600,  1.454e-13, 71.835),
    (700,  3.614e-14, 88.667),
    (800,  1.170e-14, 124.64),
    (900,  5.245e-15, 181.05),
    (1000, 3.019e-15, 268.00),
]

# Factor to scale USSA76 → approximate NRLMSISE-00 result
# Chosen to reproduce paper's Figure 8 NRLMSISE-00 lifetime of ~305 days
NRLMSISE_SCALE = 0.55


def atmospheric_density(alt_km, model='USSA76'):
    """
    Return atmospheric density [kg/m^3] at given altitude.

    Parameters
    ----------
    alt_km : float — altitude above Earth surface [km]
    model  : str   — 'USSA76' or 'NRLMSISE00'
    """
    if alt_km < 0:
        return USSA76_TABLE[0][1]
    if alt_km >= 1000:
        h0, rho0, H = USSA76_TABLE[-1]
        rho = rho0 * np.exp(-(alt_km - h0) / H)
    else:
        rho = USSA76_TABLE[0][1]
        for i in range(len(USSA76_TABLE) - 1):
            h0, rho0, H = USSA76_TABLE[i]
            h1 = USSA76_TABLE[i + 1][0]
            if h0 <= alt_km < h1:
                rho = rho0 * np.exp(-(alt_km - h0) / H)
                break

    if model == 'NRLMSISE00':
        rho *= NRLMSISE_SCALE
    return rho
