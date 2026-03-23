"""
Coordinate transformation utilities for COMTAILS simulation.

This module provides functions for transforming coordinates between
different reference frames used in the simulation.
"""
import numpy as np
from constants import FLOAT_TYPE, TORAD


def he_to_hpo(x_in, y_in, z_in, helio_matrix):
    """
    Convert heliocentric ecliptic coordinates to heliocentric plane-of-orbit.

    Args:
        x_in, y_in, z_in: Input coordinates
        helio_matrix: Transformation matrix

    Returns:
        tuple: (x_out, y_out, z_out) transformed coordinates
    """
    # Convert inputs to float64
    x_in = FLOAT_TYPE(x_in)
    y_in = FLOAT_TYPE(y_in)
    z_in = FLOAT_TYPE(z_in)

    xx, xy, xz, yx, yy, yz, zx, zy, zz = [FLOAT_TYPE(m) for m in helio_matrix]

    x_out = FLOAT_TYPE(xx * x_in + xy * y_in + xz * z_in)
    y_out = FLOAT_TYPE(yx * x_in + yy * y_in + yz * z_in)
    z_out = FLOAT_TYPE(zx * x_in + zy * y_in + zz * z_in)

    return x_out, y_out, z_out


def hpo_to_he(x_in, y_in, z_in, helio_matrix):
    """
    Convert heliocentric plane-of-orbit coordinates to heliocentric ecliptic.

    Args:
        x_in, y_in, z_in: Input coordinates
        helio_matrix: Transformation matrix

    Returns:
        tuple: (x_out, y_out, z_out) transformed coordinates
    """
    # Convert inputs to float64
    x_in = FLOAT_TYPE(x_in)
    y_in = FLOAT_TYPE(y_in)
    z_in = FLOAT_TYPE(z_in)

    xx, xy, xz, yx, yy, yz, zx, zy, zz = [FLOAT_TYPE(m) for m in helio_matrix]

    x_out = FLOAT_TYPE(xx * x_in + yx * y_in + zx * z_in)
    y_out = FLOAT_TYPE(xy * x_in + yy * y_in + zy * z_in)
    z_out = FLOAT_TYPE(xz * x_in + yz * y_in + zz * z_in)

    return x_out, y_out, z_out


def vectorial(a1, a2, a3, b1, b2, b3):
    """
    Calculate cross product C = A×B.

    Args:
        a1, a2, a3: Components of vector A
        b1, b2, b3: Components of vector B

    Returns:
        tuple: (c1, c2, c3) components of cross product
    """
    # Convert all inputs to float64
    a1, a2, a3 = FLOAT_TYPE(a1), FLOAT_TYPE(a2), FLOAT_TYPE(a3)
    b1, b2, b3 = FLOAT_TYPE(b1), FLOAT_TYPE(b2), FLOAT_TYPE(b3)

    c1 = FLOAT_TYPE(a2 * b3 - a3 * b2)
    c2 = FLOAT_TYPE(a3 * b1 - a1 * b3)
    c3 = FLOAT_TYPE(a1 * b2 - a2 * b1)

    return c1, c2, c3


def std_coor(alpha0_in, delta0_in, alpha_in, delta_in):
    """
    Convert (RA, DEC) to standard coordinates (X, Y).
    (RA, DEC) in degrees.

    Args:
        alpha0_in: Reference right ascension
        delta0_in: Reference declination
        alpha_in: Target right ascension
        delta_in: Target declination

    Returns:
        tuple: (x, y) standard coordinates
    """
    # Convert inputs to float64
    alpha0_in = FLOAT_TYPE(alpha0_in)
    delta0_in = FLOAT_TYPE(delta0_in)
    alpha_in = FLOAT_TYPE(alpha_in)
    delta_in = FLOAT_TYPE(delta_in)

    # Convert to radians
    alpha0 = FLOAT_TYPE(alpha0_in * TORAD)
    delta0 = FLOAT_TYPE(delta0_in * TORAD)
    alpha = FLOAT_TYPE(alpha_in * TORAD)
    delta = FLOAT_TYPE(delta_in * TORAD)

    deno = FLOAT_TYPE(
        np.cos(delta0) * np.cos(delta) * np.cos(alpha - alpha0) +
        np.sin(delta) * np.sin(delta0)
    )

    x = FLOAT_TYPE(np.cos(delta) * np.sin(alpha - alpha0) / deno)

    y = FLOAT_TYPE(
        (np.sin(delta0) * np.cos(delta) * np.cos(alpha - alpha0) -
         np.cos(delta0) * np.sin(delta)) / deno
    )

    return x, y


class TransformationFactory:
    """Factory class for creating coordinate transformation matrices."""

    @staticmethod
    def set_helio_matrix(inc, om, wper):
        """
        Set rotation matrix elements from heliocentric orbital plane
        to heliocentric ecliptic axes.

        Args:
            inc: Inclination (radians)
            om: Longitude of ascending node (radians)
            wper: Argument of perihelion (radians)

        Returns:
            tuple: 9 elements of the rotation matrix
        """
        # Convert inputs to float64
        inc = FLOAT_TYPE(inc)
        om = FLOAT_TYPE(om)
        wper = FLOAT_TYPE(wper)

        cw = np.cos(wper)
        co = np.cos(om)
        ci = np.cos(inc)
        sw = np.sin(wper)
        so = np.sin(om)
        si = np.sin(inc)

        xx = FLOAT_TYPE(cw * co - ci * so * sw)
        yx = FLOAT_TYPE(-sw * co - ci * so * cw)
        zx = FLOAT_TYPE(si * so)

        xy = FLOAT_TYPE(cw * so + ci * co * sw)
        yy = FLOAT_TYPE(-sw * so + ci * co * cw)
        zy = FLOAT_TYPE(-si * co)

        xz = FLOAT_TYPE(sw * si)
        yz = FLOAT_TYPE(cw * si)
        zz = FLOAT_TYPE(ci)

        return (xx, xy, xz, yx, yy, yz, zx, zy, zz)
