"""
Heliocentric orbit calculation module for COMTAILS simulation.

This module provides classes for heliocentric orbit calculations and
coordinate transformations for dust particles.
"""
import numpy as np

from constants import FLOAT_TYPE, TWOPI, MU, AUKM
from utils.coordinate_transforms import he_to_hpo
from orbital.orbit_solver import ekepl2, hkepler


class HelioOrbit:
    """
    Class for handling heliocentric orbit calculations.

    This class handles the orbital mechanics for dust particles,
    calculating their positions in various coordinate systems.
    """

    def __init__(self):
        """Initialize the heliocentric orbit handler."""
        self.helio_matrix = None
        self.helio_matrix_particle = None
        self.orbital_elements_particle = None
        self.params_nm = None

    def set_matrices(self, helio_matrix, helio_matrix_particle):
        """
        Set the transformation matrices.

        Args:
            helio_matrix: Heliocentric transformation matrix for comet
            helio_matrix_particle: Heliocentric transformation matrix for particle
        """
        self.helio_matrix = helio_matrix
        self.helio_matrix_particle = helio_matrix_particle

    def set_orbital_elements(self, qd, ed, ind, wd, nd):
        """
        Set orbital elements for the particle.

        Args:
            qd: Perihelion distance
            ed: Eccentricity
            ind: Inclination
            wd: Argument of perihelion
            nd: Longitude of ascending node
        """
        # Convert inputs to float64
        qd = FLOAT_TYPE(qd)
        ed = FLOAT_TYPE(ed)
        ind = FLOAT_TYPE(ind)
        wd = FLOAT_TYPE(wd)
        nd = FLOAT_TYPE(nd)
        self.orbital_elements_particle = (qd, ed, ind, wd, nd)

    def set_params_nm(self, delta, nmpar, nmpar1, nmpar2, nmpar3, nmpar4,
                     nmpar5, nmpar6, nmpar7, nmpar8):
        """
        Set parameters for sky plane coordinate transformation.

        Args:
            delta: Distance to observer
            nmpar, nmpar1, ...: Sky plane transformation parameters
        """
        # Convert all parameters to float64
        delta = FLOAT_TYPE(delta)
        nmpar = FLOAT_TYPE(nmpar)
        nmpar1 = FLOAT_TYPE(nmpar1)
        nmpar2 = FLOAT_TYPE(nmpar2)
        nmpar3 = FLOAT_TYPE(nmpar3)
        nmpar4 = FLOAT_TYPE(nmpar4)
        nmpar5 = FLOAT_TYPE(nmpar5)
        nmpar6 = FLOAT_TYPE(nmpar6)
        nmpar7 = FLOAT_TYPE(nmpar7)
        nmpar8 = FLOAT_TYPE(nmpar8)

        self.params_nm = (delta, nmpar, nmpar1, nmpar2, nmpar3, nmpar4,
                         nmpar5, nmpar6, nmpar7, nmpar8)

    def heliorbit(self, per_jd, tc, tau, qp, pden, rad,
                 xc_ejec, yc_ejec, zc_ejec, vxc_ejec, vyc_ejec, vzc_ejec,
                 rc_obs, thetac_obs, x, y, z, vx, vy, vz):
        """
        Calculate heliocentric orbit and sky plane coordinates.

        Args:
            per_jd: Perihelion Julian date
            tc: Observation time relative to perihelion (days)
            tau: Ejection time relative to tc (days)
            qp: Scattering efficiency for radiation pressure
            pden: Particle density (kg/m³)
            rad: Particle radius (meters)
            xc_ejec, yc_ejec, zc_ejec: Heliocentric position of comet at ejection
            vxc_ejec, vyc_ejec, vzc_ejec: Heliocentric velocity of comet at ejection
            rc_obs: Comet heliocentric distance at observation
            thetac_obs: Comet true anomaly at observation (radians)
            x, y, z: Position of particle relative to comet
            vx, vy, vz: Velocity of particle relative to comet

        Returns:
            tuple: (npar, mpar, lpar) Sky plane coordinates in km
        """
        # Convert all inputs to float64
        per_jd = FLOAT_TYPE(per_jd)
        tc = FLOAT_TYPE(tc)
        tau = FLOAT_TYPE(tau)
        qp = FLOAT_TYPE(qp)
        pden = FLOAT_TYPE(pden)
        rad = FLOAT_TYPE(rad)
        xc_ejec = FLOAT_TYPE(xc_ejec)
        yc_ejec = FLOAT_TYPE(yc_ejec)
        zc_ejec = FLOAT_TYPE(zc_ejec)
        vxc_ejec = FLOAT_TYPE(vxc_ejec)
        vyc_ejec = FLOAT_TYPE(vyc_ejec)
        vzc_ejec = FLOAT_TYPE(vzc_ejec)
        rc_obs = FLOAT_TYPE(rc_obs)
        thetac_obs = FLOAT_TYPE(thetac_obs)
        x = FLOAT_TYPE(x)
        y = FLOAT_TYPE(y)
        z = FLOAT_TYPE(z)
        vx = FLOAT_TYPE(vx)
        vy = FLOAT_TYPE(vy)
        vz = FLOAT_TYPE(vz)

        # Time of ejection
        tiem0 = FLOAT_TYPE(tc - tau + per_jd)

        # Parameters for radiation pressure
        umu = FLOAT_TYPE(1.191e-3 * qp / (2.0 * pden * rad))
        gm = FLOAT_TYPE((1.0 - umu) * MU)  # Reduced mu*GM

        # Heliocentric ecliptic position and velocity of particle
        xoutp = FLOAT_TYPE(x + xc_ejec)
        youtp = FLOAT_TYPE(y + yc_ejec)
        zoutp = FLOAT_TYPE(z + zc_ejec)

        vxoutp = FLOAT_TYPE(vx + vxc_ejec)
        vyoutp = FLOAT_TYPE(vy + vyc_ejec)
        vzoutp = FLOAT_TYPE(vz + vzc_ejec)

        # Calculate orbital elements
        if umu > 1.0:  # Repulsive hyperbola
            qd, ed, ind, nd, wd, ad, thetad, ene, t0d = \
                self.he_to_orbital_elements_r(-gm, tiem0, xoutp, youtp, zoutp,
                                            vxoutp, vyoutp, vzoutp)
        else:  # Elliptic or attractive hyperbola
            qd, ed, ind, nd, wd, ad, thetad, ene, t0d = \
                self.he_to_orbital_elements(umu, gm, tiem0, xoutp, youtp, zoutp,
                                          vxoutp, vyoutp, vzoutp)

        # Set orbital elements for later use
        self.set_orbital_elements(qd, ed, ind, wd, nd)

        # Compute sky plane coordinates
        npar, mpar, lpar = self.nm(gm, tiem0, tc, tau, rc_obs, thetac_obs, thetad, t0d, umu, ad)

        # Convert to km
        npar = FLOAT_TYPE(npar * AUKM)
        mpar = FLOAT_TYPE(mpar * AUKM)
        lpar = FLOAT_TYPE(lpar * AUKM)

        return npar, mpar, lpar

    def he_to_orbital_elements(self, umu, gm, tiem0, x, y, z, vx, vy, vz):
        """
        Calculate orbital elements from heliocentric state vectors.

        Args:
            umu: 1 - mu (radiation pressure parameter)
            gm: Gravitational parameter
            tiem0: Time of ejection
            x, y, z: Position
            vx, vy, vz: Velocity

        Returns:
            tuple: Orbital elements (q, e, i, om, w, a, anom, ene, t0d)
        """
        # Convert all inputs to float64
        umu = FLOAT_TYPE(umu)
        gm = FLOAT_TYPE(gm)
        tiem0 = FLOAT_TYPE(tiem0)
        x = FLOAT_TYPE(x)
        y = FLOAT_TYPE(y)
        z = FLOAT_TYPE(z)
        vx = FLOAT_TYPE(vx)
        vy = FLOAT_TYPE(vy)
        vz = FLOAT_TYPE(vz)

        # Angular momentum vector using vectorial product
        from utils.coordinate_transforms import vectorial
        hx, hy, hz = vectorial(x, y, z, vx, vy, vz)

        h2 = FLOAT_TYPE(hx*hx + hy*hy + hz*hz)
        h = np.sqrt(h2, dtype=FLOAT_TYPE)

        r = np.sqrt(x*x + y*y + z*z, dtype=FLOAT_TYPE)
        v2 = FLOAT_TYPE(vx*vx + vy*vy + vz*vz)
        v = np.sqrt(v2, dtype=FLOAT_TYPE)
        rv = FLOAT_TYPE(x*vx + y*vy + z*vz)

        # Semi-major axis from vis viva equation
        a = FLOAT_TYPE(1.0 / (2.0/r - v2/gm))

        # Eccentricity
        e = np.sqrt(1.0 - h2/(gm*a), dtype=FLOAT_TYPE)

        # Perihelion distance
        if e != 1.0:
            q = FLOAT_TYPE(a * (1.0 - e))
        else:
            # Parabolic case (rare)
            q = FLOAT_TYPE(h2 / (2.0 * gm))

        # Node
        om = np.arctan2(hx, -hy, dtype=FLOAT_TYPE)
        if om < 0.0:
            om = FLOAT_TYPE(om + TWOPI)

        # Inclination
        cosi = FLOAT_TYPE(hz / h)
        sini = FLOAT_TYPE((hx / h) / np.sin(om, dtype=FLOAT_TYPE))
        i = np.arctan2(sini, cosi, dtype=FLOAT_TYPE)

        # Argument of perihelion and true anomaly
        u = np.arctan2(z/sini, x*np.cos(om, dtype=FLOAT_TYPE) + y*np.sin(om, dtype=FLOAT_TYPE), dtype=FLOAT_TYPE)
        rdot = FLOAT_TYPE((x*vx + y*vy + z*vz) / r)

        # True anomaly
        anom = np.arctan2(h*rdot/gm, h2/(r*gm) - 1.0, dtype=FLOAT_TYPE)

        # Argument of perihelion
        w = FLOAT_TYPE(u - anom)
        if w < 0.0:
            w = FLOAT_TYPE(w + TWOPI)
        if w > TWOPI:
            w = FLOAT_TYPE(w % TWOPI)

        # Calculate time of perihelion passage
        if e < 1.0 and umu < 1.0:  # Elliptic trajectory
            ene = np.sqrt(gm / a**3, dtype=FLOAT_TYPE)
            temp1 = np.sqrt((1.0 + e) / (1.0 - e), dtype=FLOAT_TYPE)
            temp2 = FLOAT_TYPE(2.0 * np.arctan(
                (1.0 / temp1) * np.tan(anom / 2.0, dtype=FLOAT_TYPE), dtype=FLOAT_TYPE))  # Eccentric anomaly
            t0d = FLOAT_TYPE(tiem0 - (1.0 / ene) * (temp2 - e * np.sin(temp2, dtype=FLOAT_TYPE)))
        elif e > 1.0:  # Hyperbolic (attractive)
            temp1 = FLOAT_TYPE((e - 1.0) / (e + 1.0))
            temp2 = FLOAT_TYPE(np.sqrt(temp1, dtype=FLOAT_TYPE) * np.tan(anom / 2.0, dtype=FLOAT_TYPE))
            temp1 = FLOAT_TYPE(2.0 * np.arctanh(temp2))
            ene = np.sqrt(-gm / a**3, dtype=FLOAT_TYPE)
            t0d = FLOAT_TYPE(tiem0 - (1.0 / ene) * (e * np.sinh(temp1, dtype=FLOAT_TYPE) - temp1))
        else:  # Parabolic (very rare)
            zz = np.tan(anom / 2.0, dtype=FLOAT_TYPE)
            ene = np.sqrt(gm / (2.0 * q**3), dtype=FLOAT_TYPE)
            t0d = FLOAT_TYPE(tiem0 - (1.0 / ene) * (zz + zz**3 / 3.0))

        return q, e, i, om, w, a, anom, ene, t0d

    def he_to_orbital_elements_r(self, gm, tiem0, x, y, z, vx, vy, vz):
        """
        Calculate orbital elements for repulsive hyperbola.

        Args:
            gm: Gravitational parameter (negative for repulsive)
            tiem0: Time of ejection
            x, y, z: Position
            vx, vy, vz: Velocity

        Returns:
            tuple: Orbital elements (q, e, i, om, w, a, anom, ene, t0d)
        """
        # Convert all inputs to float64
        gm = FLOAT_TYPE(gm)
        tiem0 = FLOAT_TYPE(tiem0)
        x = FLOAT_TYPE(x)
        y = FLOAT_TYPE(y)
        z = FLOAT_TYPE(z)
        vx = FLOAT_TYPE(vx)
        vy = FLOAT_TYPE(vy)
        vz = FLOAT_TYPE(vz)

        # Angular momentum vector using vectorial product
        from utils.coordinate_transforms import vectorial
        hx, hy, hz = vectorial(x, y, z, vx, vy, vz)

        h2 = FLOAT_TYPE(hx*hx + hy*hy + hz*hz)
        h = np.sqrt(h2, dtype=FLOAT_TYPE)

        r = np.sqrt(x*x + y*y + z*z, dtype=FLOAT_TYPE)
        v2 = FLOAT_TYPE(vx*vx + vy*vy + vz*vz)
        v = np.sqrt(v2, dtype=FLOAT_TYPE)

        # Semi-major axis from vis viva equation (repulsive case)
        a = FLOAT_TYPE(1.0 / (v**2/gm + 2.0/r))
        a = FLOAT_TYPE(-abs(a))

        # Eccentricity
        e = np.sqrt(1.0 - h2/(gm*a), dtype=FLOAT_TYPE)

        # Perihelion distance (for repulsive hyperbola)
        q = FLOAT_TYPE(abs(a) * (1.0 + e))

        # Node
        om = np.arctan2(hx, -hy, dtype=FLOAT_TYPE)
        if om < 0.0:
            om = FLOAT_TYPE(om + TWOPI)

        # Inclination
        cosi = FLOAT_TYPE(hz / h)
        sini = FLOAT_TYPE((hx / h) / np.sin(om, dtype=FLOAT_TYPE))
        i = np.arctan2(sini, cosi, dtype=FLOAT_TYPE)

        # Argument of perihelion and true anomaly
        u = np.arctan2(z/sini, x*np.cos(om, dtype=FLOAT_TYPE) + y*np.sin(om, dtype=FLOAT_TYPE), dtype=FLOAT_TYPE)
        rdot = FLOAT_TYPE((x*vx + y*vy + z*vz) / r)

        # True anomaly (repulsive case)
        anom = np.arctan2(h*rdot/gm, h2/(r*gm) + 1.0, dtype=FLOAT_TYPE)

        # Argument of perihelion
        w = FLOAT_TYPE(u - anom)
        if w < 0.0:
            w = FLOAT_TYPE(w + TWOPI)
        if w > TWOPI:
            w = FLOAT_TYPE(w % TWOPI)

        # Eccentric anomaly for repulsive hyperbola
        f = FLOAT_TYPE(2.0 * np.arctanh(
            np.tan(anom/2.0, dtype=FLOAT_TYPE) * np.sqrt((e+1.0)/(e-1.0), dtype=FLOAT_TYPE)))
        ene = np.sqrt(abs(gm)/abs(a)**3, dtype=FLOAT_TYPE)

        # Time of perihelion passage
        t0d = FLOAT_TYPE(tiem0 - (e*np.sinh(f, dtype=FLOAT_TYPE) + f) / ene)

        return q, e, i, om, w, a, anom, ene, t0d

    def transform_matrix_particle(self):
        """
        Compute transformation matrix for particle orbital plane.

        Returns:
            tuple: Transformation matrix elements
        """
        qd, ed, ind, wd, nd = self.orbital_elements_particle

        cw = np.cos(wd, dtype=FLOAT_TYPE)
        sw = np.sin(wd, dtype=FLOAT_TYPE)
        ci = np.cos(ind, dtype=FLOAT_TYPE)
        si = np.sin(ind, dtype=FLOAT_TYPE)
        cn = np.cos(nd, dtype=FLOAT_TYPE)
        sn = np.sin(nd, dtype=FLOAT_TYPE)

        xx = FLOAT_TYPE(cw*cn - ci*sn*sw)
        yx = FLOAT_TYPE(-sw*cn - ci*sn*cw)
        zx = FLOAT_TYPE(si*sn)

        xy = FLOAT_TYPE(cw*sn + ci*cn*sw)
        yy = FLOAT_TYPE(-sw*sn + ci*cn*cw)
        zy = FLOAT_TYPE(-si*cn)

        xz = FLOAT_TYPE(sw*si)
        yz = FLOAT_TYPE(cw*si)
        zz = FLOAT_TYPE(ci)

        return (xx, xy, xz, yx, yy, yz, zx, zy, zz)

    def nm(self, gm, tiem0, tc_in, tau_in, rc, thetac, thetad, t0dex, umu, ad):
        """
        Calculate sky plane coordinates (N, M, L).

        Args:
            gm: Gravitational parameter
            tiem0: Time of ejection
            tc_in: Observation time
            tau_in: Ejection time relative to tc
            rc: Comet heliocentric distance
            thetac: Comet true anomaly
            thetad: Dust particle true anomaly at ejection
            t0dex: Time of dust perihelion passage
            umu: Radiation pressure parameter
            ad: Dust semi-major axis

        Returns:
            tuple: (npar, mpar, lpar) Sky plane coordinates
        """
        # Convert inputs to float64
        gm = FLOAT_TYPE(gm)
        tiem0 = FLOAT_TYPE(tiem0)
        tc = FLOAT_TYPE(tc_in)
        tau = FLOAT_TYPE(tau_in)
        rc = FLOAT_TYPE(rc)
        thetac = FLOAT_TYPE(thetac)
        thetad = FLOAT_TYPE(thetad)
        t0dex = FLOAT_TYPE(t0dex)
        umu = FLOAT_TYPE(umu)
        ad = FLOAT_TYPE(ad)

        qd, ed, ind, wd, nd = [FLOAT_TYPE(param) for param in self.orbital_elements_particle]
        delta, nmpar, nmpar1, nmpar2, nmpar3, nmpar4, nmpar5, nmpar6, nmpar7, nmpar8 = [FLOAT_TYPE(param) for param in self.params_nm]

        # Calculate position based on orbital type
        if ed < 1.0 and umu < 1.0:  # Elliptic trajectory
            ene = np.sqrt(gm / ad**3, dtype=FLOAT_TYPE)
            eme = FLOAT_TYPE(ene * (tiem0 - t0dex + tau))
            eep = FLOAT_TYPE(ekepl2(eme, ed))
            thetad = FLOAT_TYPE(2.0 * np.arctan(np.sqrt((ed+1.0)/(1.0-ed), dtype=FLOAT_TYPE) * np.tan(eep/2.0, dtype=FLOAT_TYPE), dtype=FLOAT_TYPE))
            rd = FLOAT_TYPE(ad * (1.0 - ed*ed) / (1.0 + ed*np.cos(thetad, dtype=FLOAT_TYPE)))
        else:  # Hyperbolic trajectories
            if (1.0 - umu) > 0.0:  # Attractive hyperbola
                ene = np.sqrt(-gm / ad**3, dtype=FLOAT_TYPE)
                eme = FLOAT_TYPE(ene * (tiem0 - t0dex + tau))
                eep = FLOAT_TYPE(hkepler(eme, ed))
                rd = FLOAT_TYPE(-ad * (ed*np.cosh(eep, dtype=FLOAT_TYPE) - 1.0))
                thetad = FLOAT_TYPE(2.0 * np.arctan(np.sqrt((ed+1.0)/(ed-1.0), dtype=FLOAT_TYPE) * np.tanh(eep/2.0, dtype=FLOAT_TYPE), dtype=FLOAT_TYPE))
            else:  # Repulsive hyperbola
                ene = np.sqrt(gm / ad**3, dtype=FLOAT_TYPE)
                eme = FLOAT_TYPE(ene * (tiem0 - t0dex + tau))
                eep = FLOAT_TYPE(hkepler(-eme, -ed))
                thetad = FLOAT_TYPE(2.0 * np.arctan(np.sqrt((ed-1.0)/(ed+1.0), dtype=FLOAT_TYPE) * np.tanh(eep/2.0, dtype=FLOAT_TYPE), dtype=FLOAT_TYPE))
                rd = FLOAT_TYPE(abs(ad) * (ed*np.cosh(eep, dtype=FLOAT_TYPE) + 1.0))

        # Particle position in orbital plane
        xd = FLOAT_TYPE(rd * np.cos(thetad, dtype=FLOAT_TYPE))
        yd = FLOAT_TYPE(rd * np.sin(thetad, dtype=FLOAT_TYPE))
        zd = FLOAT_TYPE(0.0)

        # Update particle transformation matrix
        self.helio_matrix_particle = self.transform_matrix_particle()

        # Convert to heliocentric ecliptic
        xde, yde, zde = self._hpo_to_he_particle(xd, yd, zd)

        # Convert to comet orbital plane
        xout, yout, zout = he_to_hpo(xde, yde, zde, self.helio_matrix)

        # Compute cometocentric coordinates
        sint = np.sin(thetac, dtype=FLOAT_TYPE)
        cost = np.cos(thetac, dtype=FLOAT_TYPE)
        chita = FLOAT_TYPE(xout*cost + yout*sint - rc)
        eta = FLOAT_TYPE(xout*sint - yout*cost)
        gita = FLOAT_TYPE(zout)

        # Convert to sky plane coordinates
        mpar = FLOAT_TYPE(nmpar1*chita - nmpar2*eta - nmpar3*gita)
        npar = FLOAT_TYPE(nmpar4*eta - nmpar5*gita)
        lpar = FLOAT_TYPE(nmpar6*chita + nmpar7*eta + nmpar8*gita)

        return npar, mpar, lpar

    def _hpo_to_he_particle(self, x_in, y_in, z_in):
        """
        Convert heliocentric plane-of-orbit to heliocentric ecliptic coordinates.

        Args:
            x_in, y_in, z_in: Input coordinates

        Returns:
            tuple: (x_out, y_out, z_out) transformed coordinates
        """
        # Convert inputs to float64
        x_in = FLOAT_TYPE(x_in)
        y_in = FLOAT_TYPE(y_in)
        z_in = FLOAT_TYPE(z_in)

        xxp, xyp, xzp, yxp, yyp, yzp, zxp, zyp, zzp = [FLOAT_TYPE(m) for m in self.helio_matrix_particle]

        x_out = FLOAT_TYPE(xxp*x_in + yxp*y_in + zxp*z_in)
        y_out = FLOAT_TYPE(xyp*x_in + yyp*y_in + zyp*z_in)
        z_out = FLOAT_TYPE(xzp*x_in + yzp*y_in + zzp*z_in)

        return x_out, y_out, z_out
