"""
Date and time utilities for COMTAILS simulation.

This module provides date/time conversion functions for astronomical calculations.
"""
from constants import FLOAT_TYPE


def caldate(jd):
    """
    Convert Julian Date to calendar date (day, month, year).
    Based on "Practical Astronomy with your calculator" by P.Duffett-Smith.

    Args:
        jd: Julian Date

    Returns:
        tuple: (day, month, year)
    """
    jd = FLOAT_TYPE(jd)
    xjd = FLOAT_TYPE(jd + 0.5)
    i = int(xjd)
    f = FLOAT_TYPE(xjd - i)

    if i > 2299160:
        ar = FLOAT_TYPE(i - 1867216.25) / FLOAT_TYPE(36524.25)
        a = int(ar)
        b = i + 1 + a - int(FLOAT_TYPE(a) / FLOAT_TYPE(4.0))
    else:
        b = i

    c = FLOAT_TYPE(b + 1524)
    d = int((FLOAT_TYPE(c) - FLOAT_TYPE(122.1)) / FLOAT_TYPE(365.25))
    e = int(FLOAT_TYPE(365.25) * FLOAT_TYPE(d))
    g = int((FLOAT_TYPE(c - e)) / FLOAT_TYPE(30.6001))

    dd = FLOAT_TYPE(c - e + f - int(FLOAT_TYPE(30.6001) * FLOAT_TYPE(g)))

    if g < 13.5:
        mm = g - 1
    else:
        mm = g - 13

    if mm > 2.5:
        yy = d - 4716
    else:
        yy = d - 4715

    return dd, mm, yy


