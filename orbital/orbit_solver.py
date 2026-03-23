"""
Orbit solver module for COMTAILS simulation.

This module provides functions for solving Kepler's equation and
computing positions and velocities from orbital elements.
"""
import numpy as np
from constants import FLOAT_TYPE, PI, TWOPI


def elements_to_xv(gm, t, qr, ec, per_jd, helio_matrix):
    """
    Convert orbital elements to position and velocity vectors.
    Enhanced for numerical stability and precision.

    Args:
        gm: Gravitational parameter
        t: Time
        qr: Perihelion distance
        ec: Eccentricity
        per_jd: Perihelion time
        helio_matrix: Heliocentric transformation matrix

    Returns:
        tuple: (x, y, z, vx, vy, vz, theta) position, velocity, and true anomaly
    """
    # Convert inputs to float64
    gm = FLOAT_TYPE(gm)
    t = FLOAT_TYPE(t)
    qr = FLOAT_TYPE(qr)
    ec = FLOAT_TYPE(ec)
    per_jd = FLOAT_TYPE(per_jd)

    if ec < 1.0:  # Elliptic trajectory
        ax = FLOAT_TYPE(qr / (1.0 - ec))
        ene = FLOAT_TYPE(np.sqrt(gm / ax ** 3))
        eme = FLOAT_TYPE(ene * (t - per_jd))

        # Use more stable method for Kepler's equation
        ee = ekepl2(eme, ec)

        # More numerically stable formula for true anomaly
        sin_e_half = FLOAT_TYPE(np.sin(ee / 2.0))
        cos_e_half = FLOAT_TYPE(np.cos(ee / 2.0))
        tan_theta_half = FLOAT_TYPE(np.sqrt((1.0 + ec) / (1.0 - ec)) * sin_e_half / cos_e_half)
        theta = FLOAT_TYPE(2.0 * np.arctan(tan_theta_half))

        # Calculate r directly from orbital parameters
        r = FLOAT_TYPE(ax * (1.0 - ec * np.cos(ee)))

    elif ec > 1.0:  # Hyperbolic trajectory
        ax = FLOAT_TYPE(qr / (ec - 1.0))  # Negative semi-major axis
        ene = FLOAT_TYPE(np.sqrt(-gm / ax ** 3))
        eme = FLOAT_TYPE(ene * (t - per_jd))

        # Solve hyperbolic Kepler equation with improved stability
        ee = hkepler(eme, ec)

        # Calculate r from the hyperbolic form
        r = FLOAT_TYPE(-ax * (ec * np.cosh(ee) - 1.0))

        # Calculate true anomaly directly from eccentric anomaly
        sinh_e_half = FLOAT_TYPE(np.sinh(ee / 2.0))
        cosh_e_half = FLOAT_TYPE(np.cosh(ee / 2.0))
        tan_theta_half = FLOAT_TYPE(np.sqrt((ec + 1.0) / (ec - 1.0)) * sinh_e_half / cosh_e_half)
        theta = FLOAT_TYPE(2.0 * np.arctan(tan_theta_half))

    else:  # Parabolic trajectory
        onethird = FLOAT_TYPE(1.0 / 3.0)
        cc = FLOAT_TYPE(3.0 * np.sqrt(gm / (2.0 * qr ** 3)) * (t - per_jd))

        # Use more stable algorithm for parabolic case
        if abs(cc) < 1.0e-10:
            zpar = FLOAT_TYPE(0.0)  # Handle very small values
        else:
            # Solve cubic equation directly
            wpar = FLOAT_TYPE(np.cbrt(np.sqrt(1.0 + (cc / 2.0) ** 2) + abs(cc / 2.0)) *
                             np.sign(cc))
            zpar = FLOAT_TYPE(wpar - 1.0 / wpar)

        theta = FLOAT_TYPE(2.0 * np.arctan(zpar))
        r = FLOAT_TYPE(qr * (1.0 + zpar ** 2))

    # Position in orbital plane
    x_orb = FLOAT_TYPE(r * np.cos(theta))
    y_orb = FLOAT_TYPE(r * np.sin(theta))
    z_orb = FLOAT_TYPE(0.0)

    # Transform to heliocentric coordinates
    from utils.coordinate_transforms import hpo_to_he
    x, y, z = hpo_to_he(x_orb, y_orb, z_orb, helio_matrix)

    # Calculate velocities
    if ec < 1.0:  # Elliptic
        vc = FLOAT_TYPE(np.sqrt(gm * ax) / r)
        # Use more precise form for velocity components
        vx_orb = FLOAT_TYPE(-vc * np.sin(ee))
        vy_orb = FLOAT_TYPE(vc * np.sqrt(1.0 - ec * ec) * np.cos(ee))
        vz_orb = FLOAT_TYPE(0.0)
    elif ec > 1.0:  # Hyperbolic
        vc = FLOAT_TYPE(-np.sqrt(-gm * ax) / r)
        h = FLOAT_TYPE(2.0 * np.arctanh(np.sqrt((ec - 1.0) / (ec + 1.0)) * np.tan(theta / 2.0)))
        vx_orb = FLOAT_TYPE(vc * np.sinh(h))
        vy_orb = FLOAT_TYPE(-vc * np.sqrt(ec * ec - 1.0) * np.cosh(h))
        vz_orb = FLOAT_TYPE(0.0)
    else:  # Parabolic
        # More stable formula for parabolic velocity
        vy_orb = FLOAT_TYPE(np.sqrt(2.0 * gm / qr) / (1.0 + zpar ** 2))
        vx_orb = FLOAT_TYPE(-zpar * vy_orb)
        vz_orb = FLOAT_TYPE(0.0)

    # Transform velocities to heliocentric
    vx, vy, vz = hpo_to_he(vx_orb, vy_orb, vz_orb, helio_matrix)

    return x, y, z, vx, vy, vz, theta


def ekepl2(em, e):
    """
    Solves Kepler equation EM = E - e*sin(E) for elliptical orbits.

    Args:
        em: Mean anomaly
        e: Eccentricity

    Returns:
        float: Eccentric anomaly
    """
    # Convert inputs to float64
    em = FLOAT_TYPE(em)
    e = FLOAT_TYPE(e)

    # Port of the original EKEPL2 function
    pineg = FLOAT_TYPE(-PI)
    sw = FLOAT_TYPE(0.1)
    ahalf = FLOAT_TYPE(0.5)
    asixth = FLOAT_TYPE(ahalf / 3.0)
    athird = FLOAT_TYPE(asixth * 2.0)
    a = FLOAT_TYPE((PI - 1.0)**2 / (PI + 2.0/3.0))
    b = FLOAT_TYPE(2.0 * (PI - asixth)**2 / (PI + 2.0/3.0))

    emr = FLOAT_TYPE(em % TWOPI)
    if emr < pineg:
        emr = FLOAT_TYPE(emr + TWOPI)
    if emr > PI:
        emr = FLOAT_TYPE(emr - TWOPI)

    ee = FLOAT_TYPE(emr)

    if ee < 0.0:
        ee = FLOAT_TYPE(-ee)
        neg_flag = True
    else:
        neg_flag = False

    if ee < asixth:
        ee = FLOAT_TYPE((6.0 * ee)**athird)
    else:
        w = FLOAT_TYPE(PI - ee)
        ee = FLOAT_TYPE(PI - a*w/(b-w))

    if neg_flag:
        ee = FLOAT_TYPE(-ee)

    ee = FLOAT_TYPE(emr + (ee - emr) * e)
    e1 = FLOAT_TYPE(1.0 - e)
    l = (e1 + ee*ee/6.0) >= sw

    for _ in range(2):
        fdd = FLOAT_TYPE(e * np.sin(ee))
        fddd = FLOAT_TYPE(e * np.cos(ee))

        if l:
            f = FLOAT_TYPE((ee - fdd) - emr)
            fd = FLOAT_TYPE(1.0 - fddd)
        else:
            f = FLOAT_TYPE(em_kepl(e, ee) - emr)
            fd = FLOAT_TYPE(e1 + 2.0 * e * np.sin(ahalf * ee)**2)

        dee = FLOAT_TYPE(f * fd / (0.5 * f * fd - fd * fd))
        w = FLOAT_TYPE(fd + 0.5 * dee * (fdd + athird * dee * fddd))
        fd = FLOAT_TYPE(fd + dee * (fdd + 0.5 * dee * fddd))
        ee = FLOAT_TYPE(ee - (f - dee * (fd - w)) / fd)

    return FLOAT_TYPE(ee + (em - emr))


def em_kepl(e, ee):
    """
    Auxiliary function for Kepler equation solver.

    Args:
        e: Eccentricity
        ee: Eccentric anomaly

    Returns:
        float: Value for mean anomaly calculation
    """
    # Convert inputs to float64
    e = FLOAT_TYPE(e)
    ee = FLOAT_TYPE(ee)

    x = FLOAT_TYPE((1.0 - e) * np.sin(ee))
    ee2 = FLOAT_TYPE(-ee * ee)
    term = FLOAT_TYPE(ee)
    d = FLOAT_TYPE(0.0)

    while True:
        d = FLOAT_TYPE(d + 2.0)
        term = FLOAT_TYPE(term * ee2 / (d * (d + 1.0)))
        x0 = x
        x = FLOAT_TYPE(x - term)
        if x == x0:
            break

    return x


def hkepler(eme, e):
    """
    Solves hyperbolic Kepler equation M = e*sinh(H)-H.

    Args:
        eme: Mean anomaly
        e: Eccentricity

    Returns:
        float: Hyperbolic anomaly
    """
    # Convert inputs to float64
    eme = FLOAT_TYPE(eme)
    e = FLOAT_TYPE(e)

    el = FLOAT_TYPE(eme / e)
    g1 = FLOAT_TYPE(1.0 - 1.0 / e)
    return np.arcsinh(sh_kepl(el, g1))


def sh_kepl(el, g1):
    """
    Solves equation: EL = SHKEPL + (G1-1) * ASINH(SHKEPL).

    Args:
        el: Parameter
        g1: Parameter

    Returns:
        float: Solution
    """
    # Convert inputs to float64
    el = FLOAT_TYPE(el)
    g1 = FLOAT_TYPE(g1)

    sw = FLOAT_TYPE(0.5)
    ahalf = FLOAT_TYPE(0.5)
    asixth = FLOAT_TYPE(ahalf / 3.0)
    athird = FLOAT_TYPE(asixth * 2.0)

    if el == 0.0:
        return FLOAT_TYPE(0.0)

    g = FLOAT_TYPE(1.0 - g1)
    cl = np.sqrt(1.0 + el**2)
    al = np.arcsinh(el)

    w = FLOAT_TYPE(g**2 * al / cl**3)
    s = FLOAT_TYPE(1.0 - g / cl)
    s = FLOAT_TYPE(el + g * al / np.cbrt(s**3 + w * el * (1.5 - g / 0.75)))

    # Two iterations of halley-then-newton process
    for _ in range(2):
        s0 = FLOAT_TYPE(s * s)
        s1 = FLOAT_TYPE(s0 + 1.0)
        s2 = np.sqrt(s1)
        s3 = FLOAT_TYPE(s1 * s2)
        fdd = FLOAT_TYPE(g * s / s3)
        fddd = FLOAT_TYPE(g * (1.0 - 2.0 * s0) / (s1 * s3))

        if asixth * s0 + g1 >= sw:
            f = FLOAT_TYPE((s - g * np.arcsinh(s)) - el)
            fd = FLOAT_TYPE(1.0 - g / s2)
        else:
            f = FLOAT_TYPE(shm_kep(g1, s) - el)
            fd = FLOAT_TYPE((s0 / (s2 + 1.0) + g1) / s2)

        ds = FLOAT_TYPE(f * fd / (ahalf * f * fdd - fd * fd))
        stemp = FLOAT_TYPE(s + ds)

        if stemp == s:
            break

        f = FLOAT_TYPE(f + ds * (fd + ahalf * ds * (fdd + athird * ds * fddd)))
        fd = FLOAT_TYPE(fd + ds * (fdd + ahalf * ds * fddd))
        s = FLOAT_TYPE(stemp - f / fd)

    return s


def shm_kep(g1, s):
    """
    Helper function for hyperbolic Kepler equation.

    Args:
        g1: Parameter
        s: Parameter

    Returns:
        float: Result
    """
    # Convert inputs to float64
    g1 = FLOAT_TYPE(g1)
    s = FLOAT_TYPE(s)

    g = FLOAT_TYPE(1.0 - g1)
    t = FLOAT_TYPE(s / (1.0 + np.sqrt(1.0 + s * s)))
    tsq = FLOAT_TYPE(t * t)
    x = FLOAT_TYPE(s * (g1 + g * tsq))
    term = FLOAT_TYPE(2.0 * g * t)
    twoi1 = FLOAT_TYPE(1.0)

    while True:
        twoi1 = FLOAT_TYPE(twoi1 + 2.0)
        term = FLOAT_TYPE(term * tsq)
        x0 = x
        x = FLOAT_TYPE(x - term / twoi1)
        if x == x0:
            break

    return x
