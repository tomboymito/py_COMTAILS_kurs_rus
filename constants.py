"""
Global constants for the COMTAILS simulation.

This module defines all physical and mathematical constants needed for the simulation,
ensuring consistent values throughout the codebase.
"""
import numpy as np

# Define high-precision floating-point type
# Use float64 for all calculations
FLOAT_TYPE = np.float64

# Mathematical constants
PI = FLOAT_TYPE(2.0 * np.arccos(0.0, dtype=FLOAT_TYPE))
HALFPI = FLOAT_TYPE(PI / 2.0)
TWOPI = FLOAT_TYPE(2.0 * PI)
FOURPI = FLOAT_TYPE(4.0 * PI)
TORAD = FLOAT_TYPE(PI / 180.0)

# Astronomical constants
AUKM = FLOAT_TYPE(1.4959787e8)        # Astronomical unit (km)
CTEVEL = FLOAT_TYPE(5.775483e-4)      # km/s to AU/day conversion factor
MU = FLOAT_TYPE(2.959122082855911e-4) # Standard gravitational parameter GM au³/day²
QP = FLOAT_TYPE(1.0)                  # Scattering efficiency for radiation pressure
RSUNKM = FLOAT_TYPE(695660.0)         # Solar radius, km

# Schleicher phase function polynomial fit coefficients
SCHLEICHER_COEFFS = np.array([
    -7.4004978e-03, -1.6492566e-02, 1.0950353e-04,
    8.3640600e-07, 1.0157539e-09, -9.6882641e-11, 4.4184372e-13
], dtype=FLOAT_TYPE)


