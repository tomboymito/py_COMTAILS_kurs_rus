"""
Numerical utilities for COMTAILS simulation.

This module provides numerical functions for high-precision calculations.
"""
import numpy as np
from constants import FLOAT_TYPE


def powerint(xa, xb, alpha):
    """
    Integral of a single power-law distribution of index alpha
    between radii xa and xb.

    Args:
        xa: Lower bound
        xb: Upper bound
        alpha: Power-law index

    Returns:
        float: Integral value
    """
    # Convert inputs to float64
    xa = FLOAT_TYPE(xa)
    xb = FLOAT_TYPE(xb)
    alpha = FLOAT_TYPE(alpha)

    if alpha != -1.0:
        ualp = FLOAT_TYPE(1.0 + alpha)
        return FLOAT_TYPE(xb**ualp - xa**ualp) / ualp
    else:
        return FLOAT_TYPE(np.log(xb) - np.log(xa))


