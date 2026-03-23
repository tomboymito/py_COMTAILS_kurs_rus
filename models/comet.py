"""
Comet model for COMTAILS simulation.

This module provides a class for representing a comet and its orbital properties.
"""
import numpy as np
from constants import FLOAT_TYPE, TORAD, TWOPI, MU
from utils.date_time import caldate
from orbital.orbit_solver import elements_to_xv


class Comet:
    """
    Class representing a comet and its orbital properties.

    This class stores the comet's physical parameters, orbital elements,
    and state vectors over time.
    """

    def __init__(self, config):
        """
        Initialize comet object.

        Args:
            config: SimulationConfig object containing configuration parameters
        """
        # Basic identification
        self.recid = config.recid

        # Physical parameters
        self.rnucleus = config.rnucleus
        self.pv0_nuc = config.pv0_nuc
        self.phase_coeff_nuc = config.phase_coeff_nuc
        self.brnucleus = FLOAT_TYPE(0.0)

        # Rotation parameters
        self.nuc_inc = config.nuc_inc
        self.nuc_phi = config.nuc_phi
        self.period = config.period

        # Orbital elements - will be set later
        self.ec = FLOAT_TYPE(0.5)  # eccentricity
        self.qr = FLOAT_TYPE(1.0)  # perihelion distance (AU)
        self.om = FLOAT_TYPE(0.0)  # longitude of ascending node (deg)
        self.wper = FLOAT_TYPE(0.0)  # argument of perihelion (deg)
        self.inc = FLOAT_TYPE(0.0)  # inclination (deg)

        self.per_jd = FLOAT_TYPE(0.0)

        # Orbit parameters in radians
        self.om_rad = FLOAT_TYPE(0.0)
        self.wper_rad = FLOAT_TYPE(0.0)
        self.inc_rad = FLOAT_TYPE(0.0)

        # Derived orbital parameters
        self.semiaxis = FLOAT_TYPE(0.0)
        self.orb_per = FLOAT_TYPE(0.0)
        self.half_per = FLOAT_TYPE(0.0)

        # State vectors
        self.xarr = None
        self.yarr = None
        self.zarr = None
        self.vxarr = None
        self.vyarr = None
        self.vzarr = None
        self.trueanarr = None

        # Observational parameters
        self.ra = FLOAT_TYPE(0.0)  # Right ascension (deg)
        self.dec = FLOAT_TYPE(0.0)  # Declination (deg)
        self.delta = FLOAT_TYPE(1.0)  # Distance to observer (AU)
        self.deldot = FLOAT_TYPE(0.0)  # Rate of change of distance
        self.psang = FLOAT_TYPE(0.0)  # Position angle
        self.psamv = FLOAT_TYPE(0.0)  # Position angle of heliocentric velocity
        self.plang = FLOAT_TYPE(0.0)  # Phase angle of orbital plane
        self.phas_ang = FLOAT_TYPE(0.0)  # Phase angle

    def set_orbital_elements(self, comet_data):
        """
        Set orbital elements from JPL Horizons data.

        Args:
            comet_data: Dictionary containing orbital elements and other parameters
        """
        # Orbital elements
        self.ec = FLOAT_TYPE(comet_data.get('ec', 0.5))
        self.qr = FLOAT_TYPE(comet_data.get('qr', 1.0))
        self.om = FLOAT_TYPE(comet_data.get('om', 0.0))
        self.wper = FLOAT_TYPE(comet_data.get('wper', 0.0))
        self.inc = FLOAT_TYPE(comet_data.get('inc', 0.0))

        self.per_jd = FLOAT_TYPE(comet_data.get('per_jd', 0.0))

        # Convert to radians
        self.om_rad = FLOAT_TYPE(self.om * TORAD)
        self.wper_rad = FLOAT_TYPE(self.wper * TORAD)
        self.inc_rad = FLOAT_TYPE(self.inc * TORAD)

        # Astrometric data
        self.ra = FLOAT_TYPE(comet_data.get('ra', 0.0))
        self.dec = FLOAT_TYPE(comet_data.get('dec', 0.0))
        self.delta = FLOAT_TYPE(comet_data.get('delta', 1.0))
        self.deldot = FLOAT_TYPE(comet_data.get('deldot', 0.0))
        self.psang = FLOAT_TYPE(comet_data.get('psang', 0.0))
        self.psamv = FLOAT_TYPE(comet_data.get('psamv', 0.0))
        self.plang = FLOAT_TYPE(comet_data.get('plang', 0.0))
        self.phas_ang = FLOAT_TYPE(comet_data.get('phas_ang', 0.0))

        # Calculate derived orbital parameters
        if self.ec < 1.0:  # Elliptical orbit
            self.semiaxis = FLOAT_TYPE(self.qr / (1.0 - self.ec))
            self.orb_per = FLOAT_TYPE(TWOPI * np.sqrt(self.semiaxis ** 3 / MU))
            self.half_per = FLOAT_TYPE(0.5 * self.orb_per)
        else:  # Parabolic or hyperbolic orbit
            self.semiaxis = FLOAT_TYPE(self.qr / (self.ec - 1.0) if self.ec > 1.0 else 0.0)
            self.orb_per = FLOAT_TYPE(0.0)
            self.half_per = FLOAT_TYPE(0.0)

    def compute_positions_and_velocities(self, times, helio_matrix):
        """
        Compute comet positions and velocities for all times.

        Args:
            times: Array of times (JD)
            helio_matrix: Heliocentric transformation matrix
        """
        n_times = len(times)

        # Initialize state vector arrays
        self.xarr = np.zeros(n_times, dtype=FLOAT_TYPE)
        self.yarr = np.zeros(n_times, dtype=FLOAT_TYPE)
        self.zarr = np.zeros(n_times, dtype=FLOAT_TYPE)
        self.vxarr = np.zeros(n_times, dtype=FLOAT_TYPE)
        self.vyarr = np.zeros(n_times, dtype=FLOAT_TYPE)
        self.vzarr = np.zeros(n_times, dtype=FLOAT_TYPE)
        self.trueanarr = np.zeros(n_times, dtype=FLOAT_TYPE)

        # Compute positions and velocities for all times
        for i in range(n_times):
            x, y, z, vx, vy, vz, trueanomaly = elements_to_xv(
                MU, times[i], self.qr, self.ec, self.per_jd, helio_matrix)

            self.trueanarr[i] = FLOAT_TYPE(trueanomaly)
            self.xarr[i] = FLOAT_TYPE(x)
            self.yarr[i] = FLOAT_TYPE(y)
            self.zarr[i] = FLOAT_TYPE(z)
            self.vxarr[i] = FLOAT_TYPE(vx)
            self.vyarr[i] = FLOAT_TYPE(vy)
            self.vzarr[i] = FLOAT_TYPE(vz)

    def print_orbit_info(self, start_jd, end_jd, per_jd):
        """
        Print orbit information.

        Args:
            start_jd: Start Julian date
            end_jd: End Julian date
            per_jd: Perihelion Julian date
        """
        start_dd, start_mm, start_yy = caldate(start_jd)
        obs_dd, obs_mm, obs_yy = caldate(end_jd)

        print(f" Start activity date (Day,Month,Year): {start_dd:7.3f} {int(start_mm):2d} {int(start_yy):5d}")
        print(f"  Observing date      (Day,Month,Year): {obs_dd:7.3f} {int(obs_mm):2d} {int(obs_yy):5d}")

        print(f" JD at activation: {start_jd:14.5f}  ,i.e., {start_jd - per_jd:10.2f} days to perihelion")
        print(f" JD at observation: {end_jd:14.5f}  ,i.e., {end_jd - per_jd:10.2f} days to perihelion")

        if self.ec < 1.0:
            print(f" Orbital period: {self.orb_per / 365.25:10.3f} years")
