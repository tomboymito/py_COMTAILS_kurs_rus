"""
Configuration management for COMTAILS simulation.

This module handles all the configuration settings and parameters required
for the comet dust tail simulation.
"""
import os
import numpy as np

from constants import (
    PI, TORAD, AUKM, RSUNKM, FLOAT_TYPE,
    SCHLEICHER_COEFFS
)
from utils.coordinate_transforms import he_to_hpo
from visualization.plot_handler import PlotHandler


class SimulationConfig:
    """
    Class for managing simulation configuration parameters.

    This class handles reading input files, storing parameters,
    and managing the simulation state.
    """

    def __init__(self):
        """Initialize configuration with default values."""
        # Identification
        self.recid = ""

        # Time parameters
        self.start_jd = FLOAT_TYPE(0.0)
        self.end_jd = FLOAT_TYPE(0.0)
        self.per_jd = FLOAT_TYPE(0.0)
        self.tc = FLOAT_TYPE(0.0)

        self.times = None
        self.deltat = FLOAT_TYPE(0.0)

        # Grid parameters
        self.nx = 0
        self.ny = 0
        self.scale = FLOAT_TYPE(0.0)
        self.grdsiz = FLOAT_TYPE(0.0)
        self.inuc = 0
        self.jnuc = 0
        self.nmin = FLOAT_TYPE(0.0)
        self.nmax = FLOAT_TYPE(0.0)
        self.mmin = FLOAT_TYPE(0.0)
        self.mmax = FLOAT_TYPE(0.0)
        self.ngrid = None
        self.mgrid = None
        self.augx = FLOAT_TYPE(0.0)
        self.augy = FLOAT_TYPE(0.0)
        self.areapix = FLOAT_TYPE(0.0)

        # Particle parameters
        self.pden = FLOAT_TYPE(0.0)
        self.pv0 = FLOAT_TYPE(0.0)
        self.phase_coeff = FLOAT_TYPE(0.0)
        self.v0 = FLOAT_TYPE(0.0)
        self.gamma = FLOAT_TYPE(0.0)
        self.kappa = FLOAT_TYPE(0.0)
        self.expocos = FLOAT_TYPE(0.0)

        # Nucleus parameters
        self.nuc_inc = FLOAT_TYPE(0.0)
        self.nuc_phi = FLOAT_TYPE(0.0)
        self.period = FLOAT_TYPE(0.0)
        self.lat_min = FLOAT_TYPE(0.0)
        self.lat_max = FLOAT_TYPE(0.0)
        self.lon_min = FLOAT_TYPE(0.0)
        self.lon_max = FLOAT_TYPE(0.0)
        self.lat_min_deg = FLOAT_TYPE(0.0)
        self.lat_max_deg = FLOAT_TYPE(0.0)
        self.lon_min_deg = FLOAT_TYPE(0.0)
        self.lon_max_deg = FLOAT_TYPE(0.0)
        self.rnucleus = FLOAT_TYPE(0.0)
        self.pv0_nuc = FLOAT_TYPE(0.0)
        self.phase_coeff_nuc = FLOAT_TYPE(0.0)

        # Simulation control parameters
        self.iejec_mode = 0
        self.isun = 0
        self.ntimes = 0
        self.nsizes = 0
        self.nevent = 0
        self.ninputs = 0
        self.ntotmc = 0
        self.iconv = 0
        self.sfwhm = FLOAT_TYPE(0.0)
        self.istar = 0
        self.maglim = FLOAT_TYPE(0.0)
        self.iprn = 0
        self.igrapho = 0
        self.pcp = FLOAT_TYPE(0.0)

        # Aperture parameters
        self.iap = 0
        self.rho_ap = FLOAT_TYPE(0.0)
        self.rho_ap_orig = FLOAT_TYPE(0.0)

          # Output parameters
        self.cte_mag = FLOAT_TYPE(0.0)
        self.magsun = FLOAT_TYPE(0.0)

        # Dust loss rate parameters
        self.dtime = []
        self.dmdtlog = []
        self.velfac = []
        self.powera = []
        self.radiomin = []
        self.radiomax = []

        # Position and velocity arrays
        self.xyz_earth = [FLOAT_TYPE(0.0), FLOAT_TYPE(0.0), FLOAT_TYPE(0.0)]
        self.xe = FLOAT_TYPE(0.0)
        self.ye = FLOAT_TYPE(0.0)
        self.ze = FLOAT_TYPE(0.0)
        self.xcobs = FLOAT_TYPE(0.0)
        self.ycobs = FLOAT_TYPE(0.0)
        self.zcobs = FLOAT_TYPE(0.0)
        self.thetacobs = FLOAT_TYPE(0.0)
        self.rcobs = FLOAT_TYPE(0.0)

        # Observational parameters
        self.ra = FLOAT_TYPE(0.0)
        self.dec = FLOAT_TYPE(0.0)
        self.delta = FLOAT_TYPE(0.0)
        self.deldot = FLOAT_TYPE(0.0)
        self.psang = FLOAT_TYPE(0.0)
        self.psamv = FLOAT_TYPE(0.0)
        self.plang = FLOAT_TYPE(0.0)
        self.phas_ang = FLOAT_TYPE(0.0)

        self.rah = 0
        self.ram = FLOAT_TYPE(0.0)
        self.ras = FLOAT_TYPE(0.0)
        self.ded = 0
        self.dem = FLOAT_TYPE(0.0)
        self.des = FLOAT_TYPE(0.0)


        # Coordinate system parameters
        self.arcsec_rad = FLOAT_TYPE(180.0 * 3600.0 / PI)
        self.solid_angle = FLOAT_TYPE(PI * (RSUNKM/AUKM)**2)
        self.helio_matrix = None

        # Cometocentric parameters
        self.xoute = FLOAT_TYPE(0.0)
        self.youte = FLOAT_TYPE(0.0)
        self.zoute = FLOAT_TYPE(0.0)
        self.chitae = FLOAT_TYPE(0.0)
        self.etae = FLOAT_TYPE(0.0)
        self.gitae = FLOAT_TYPE(0.0)
        self.nmpar = FLOAT_TYPE(0.0)
        self.nmpar1 = FLOAT_TYPE(0.0)
        self.nmpar2 = FLOAT_TYPE(0.0)
        self.nmpar3 = FLOAT_TYPE(0.0)
        self.nmpar4 = FLOAT_TYPE(0.0)
        self.nmpar5 = FLOAT_TYPE(0.0)
        self.nmpar6 = FLOAT_TYPE(0.0)
        self.nmpar7 = FLOAT_TYPE(0.0)
        self.nmpar8 = FLOAT_TYPE(0.0)

        # Particle size parameters
        self.r1log = FLOAT_TYPE(0.0)
        self.r2log = FLOAT_TYPE(0.0)
        self.steprlog = FLOAT_TYPE(0.0)

        # Output arrays
        self.flux = None
        self.opt_depth = None
        self.flux_star = None

        # Output file handles
        self.dustloss_file = None
        self.star_file = None
        self.nm_file = None
        self.subsolar_file = None

        if self.istar == 1:
            self.star_file = open('output/starpos.dat', 'w')
        if self.iprn == 1:
            self.nm_file = open('output/nm.dat', 'w')
        if self.iejec_mode == 3:
            self.subsolar_file = open('output/subsolar.dat', 'w')
        self.dustloss_file = open('output/dustlossrate.dat', 'w')

    def __getstate__(self):
        # Copy the object's state (all attributes)
        state = self.__dict__.copy()
        # Set file handle attributes to None to exclude them from pickling
        state['star_file'] = None
        state['nm_file'] = None
        state['subsolar_file'] = None
        state['dustloss_file'] = None
        return state

    def __setstate__(self, state):
        # Restore the state from the pickled data
        self.__dict__.update(state)
        # Ensure file handles are None in the worker process
        self.star_file = None
        self.nm_file = None
        self.subsolar_file = None
        self.dustloss_file = None

    def read_inputs(self, required_files):
        """
        Read input parameters from configuration files.

        Args:
            required_files: List of input file paths
        """
        # Helper function to parse input lines with comments
        def parse_line(line):
            # Split at the comment marker '!'
            if '!' in line:
                value_part = line.split('!')[0].strip()
            else:
                value_part = line.strip()
            return value_part

        # Read dust loss rate parameters
        with open(required_files[1], 'r') as f:
            next(f)  # Skip blank line
            for line in f:
                parts = line.split()
                if len(parts) >= 6:
                    self.dtime.append(FLOAT_TYPE(float(parts[0])))
                    self.dmdtlog.append(FLOAT_TYPE(float(parts[1])))
                    self.velfac.append(FLOAT_TYPE(float(parts[2])))
                    self.powera.append(FLOAT_TYPE(float(parts[3])))
                    self.radiomin.append(FLOAT_TYPE(float(parts[4])))
                    self.radiomax.append(FLOAT_TYPE(float(parts[5])))

        self.ninputs = len(self.dtime)

        # Read main input file
        with open(required_files[0], 'r') as f:
            self.recid = parse_line(f.readline())
            self.pden = FLOAT_TYPE(float(parse_line(f.readline())))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.pv0 = FLOAT_TYPE(float(parts[0]))
            self.phase_coeff = FLOAT_TYPE(float(parts[1]))

            self.iejec_mode = int(parse_line(f.readline()))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.nuc_inc = FLOAT_TYPE(float(parts[0]) * TORAD)
            self.nuc_phi = FLOAT_TYPE(float(parts[1]) * TORAD)
            self.period = FLOAT_TYPE(float(parts[2]))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.lat_min_deg = FLOAT_TYPE(float(parts[0]))
            self.lat_max_deg = FLOAT_TYPE(float(parts[1]))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.lon_min_deg = FLOAT_TYPE(float(parts[0]))
            self.lon_max_deg = FLOAT_TYPE(float(parts[1]))

            # Convert to radians
            self.lon_min = FLOAT_TYPE(self.lon_min_deg * TORAD)
            self.lon_max = FLOAT_TYPE(self.lon_max_deg * TORAD)
            self.lat_min = FLOAT_TYPE(self.lat_min_deg * TORAD)
            self.lat_max = FLOAT_TYPE(self.lat_max_deg * TORAD)

            self.isun = int(parse_line(f.readline()))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.v0 = FLOAT_TYPE(float(parts[0]))
            self.gamma = FLOAT_TYPE(float(parts[1]))
            self.kappa = FLOAT_TYPE(float(parts[2]))

            self.expocos = FLOAT_TYPE(float(parse_line(f.readline())))
            self.magsun = FLOAT_TYPE(float(parse_line(f.readline())))
            self.start_jd = FLOAT_TYPE(float(parse_line(f.readline())))
            self.end_jd = FLOAT_TYPE(float(parse_line(f.readline())))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.nx = int(parts[0])
            self.ny = int(parts[1])

            self.scale = FLOAT_TYPE(float(parse_line(f.readline())))
            self.ntimes = int(parse_line(f.readline()))
            self.nsizes = int(parse_line(f.readline()))
            self.nevent = int(parse_line(f.readline()))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.iap = int(parts[0])
            self.rho_ap = FLOAT_TYPE(float(parts[1]))

            # Store original aperture value
            self.rho_ap_orig = self.rho_ap

            self.rnucleus = FLOAT_TYPE(float(parse_line(f.readline())))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.pv0_nuc = FLOAT_TYPE(float(parts[0]))
            self.phase_coeff_nuc = FLOAT_TYPE(float(parts[1]))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.iconv = int(parts[0])
            self.sfwhm = FLOAT_TYPE(float(parts[1]))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.istar = int(parts[0])
            self.maglim = FLOAT_TYPE(float(parts[1]))

            self.iprn = int(parse_line(f.readline()))

            line = parse_line(f.readline())
            parts = line.split(',')
            self.igrapho = int(parts[0])
            self.pcp = FLOAT_TYPE(float(parts[1]))

            self.pcp = FLOAT_TYPE(1.0 - self.pcp / 100.0)

            # Initialize arrays with explicit data types
            self.flux = np.zeros((self.nx, self.ny), dtype=FLOAT_TYPE)
            self.opt_depth = np.zeros((self.nx, self.ny), dtype=FLOAT_TYPE)
            self.flux_star = np.zeros((self.nx, self.ny), dtype=FLOAT_TYPE)
            self.ngrid = np.zeros(self.nx, dtype=FLOAT_TYPE)
            self.mgrid = np.zeros(self.ny, dtype=FLOAT_TYPE)

            # Initialize other parameters
            self.ntotmc = self.nsizes * self.ntimes * self.nevent

            # Prepare output files
            self._prepare_output_files()

            # Initialize astronomical constants
            self._initialize_constants()

    def _prepare_output_files(self):
        """Prepare output files for the simulation."""
        os.makedirs('output', exist_ok=True)

        if self.istar == 1:
            self.star_file = open('output/starpos.dat', 'w')
        if self.iprn == 1:
            self.nm_file = open('output/nm.dat', 'w')
        if self.iejec_mode == 3:
            self.subsolar_file = open('output/subsolar.dat', 'w')
        self.dustloss_file = open('output/dustlossrate.dat', 'w')

    def _initialize_constants(self):
        """Initialize astronomical and physical constants."""
        # Solid angle of sun disk from Earth, sr
        self.solid_angle = FLOAT_TYPE(PI * (RSUNKM/AUKM)**2)

        # arcsec per radian
        self.arcsec_rad = FLOAT_TYPE(180.0 * 3600.0 / PI)

        # solid angle in arcsec**2
        self.solid_angle = FLOAT_TYPE(self.solid_angle * self.arcsec_rad**2)

        # Magnitude constant
        self.cte_mag = FLOAT_TYPE(2.5 * np.log10(self.solid_angle))

    def setup_time_array(self):
        """Set up time array for simulation with high precision."""
        self.times = np.zeros(self.ntimes, dtype=FLOAT_TYPE)
        self.deltat = FLOAT_TYPE((self.end_jd - self.start_jd) / float(self.ntimes - 1))

        self.times[0] = self.start_jd
        for i in range(1, self.ntimes):
            self.times[i] = FLOAT_TYPE(self.times[i-1] + self.deltat)

        # Ensure exact value for the end time to avoid floating point drift
        self.times[self.ntimes-1] = self.end_jd

    def setup_image_grid(self, comet):
        """Set up dust tail image grid parameters with high precision."""

        # Copy comet info
        self.ra = comet.ra
        self.dec = comet.dec
        self.delta = comet.delta
        self.deldot = comet.deldot
        self.psang = comet.psang
        self.psamv = comet.psamv
        self.plang = comet.plang
        self.phas_ang = comet.phas_ang

        self.per_jd = comet.per_jd
        self.tc = FLOAT_TYPE(self.end_jd - self.per_jd)

        # Calculate grid parameters
        self.grdsiz = FLOAT_TYPE(self.scale * self.delta * AUKM / self.arcsec_rad)  # km/px
        gr2 = FLOAT_TYPE(self.grdsiz / 2.0)

        self.inuc = self.nx // 2
        self.jnuc = self.ny // 2

        self.nmin = FLOAT_TYPE(-((self.inuc-1) * self.grdsiz + gr2))
        self.nmax = FLOAT_TYPE((self.nx-self.inuc) * self.grdsiz - gr2)

        self.mmin = FLOAT_TYPE(-((self.jnuc-1) * self.grdsiz + gr2))
        self.mmax = FLOAT_TYPE((self.ny-self.jnuc) * self.grdsiz - gr2)

        # Setup grid arrays - ensure consistent precision
        self.ngrid[0] = FLOAT_TYPE(self.nmin + gr2)
        self.mgrid[0] = FLOAT_TYPE(self.mmin + gr2)

        # Fix the aperture radius calculation
        if self.iap == 1:
            # Aperture in arcsec needs to be converted to km
            self.rho_ap = FLOAT_TYPE((self.rho_ap_orig / self.scale) * self.grdsiz)  # Aperture in km
        else:
            self.rho_ap = FLOAT_TYPE(self.rho_ap)  # Already in km

        for i in range(1, self.nx):
            self.ngrid[i] = FLOAT_TYPE(self.ngrid[i-1] + self.grdsiz)

        for i in range(1, self.ny):
            self.mgrid[i] = FLOAT_TYPE(self.mgrid[i-1] + self.grdsiz)

        # Calculate additional parameters
        self.augx = FLOAT_TYPE(float(self.nx) / (self.nmax - self.nmin))
        self.augy = FLOAT_TYPE(float(self.ny) / (self.mmax - self.mmin))
        self.areapix = FLOAT_TYPE((self.grdsiz * self.arcsec_rad / (self.delta * AUKM))**2)

        # Initialize plot handler if plotting is enabled
        plot_handler = None
        if self.igrapho == 1:
            plot_filename = 'output/dust_particles.png'
            plot_handler = PlotHandler(
                self.nmin, self.nmax,
                self.mmin, self.mmax,
                plot_filename
            )
            # Calculate how many particles we expect to plot
            expected_count = int(self.ntotmc * (1.0 - self.pcp))
            print(f"Will plot approximately {expected_count} particles ({(1.0 - self.pcp) * 100:.1f}% of total)")

        # Print grid information
        print(f" Image scale={self.grdsiz:10.3f} km/px <==> {self.scale:8.3f} arcsec/px " +
              f"Field of view= {self.nx*self.scale/60.0:8.3f} x{self.ny*self.scale/60.0:8.3f} arcmin")

        return plot_handler

    def set_earth_position(self, position_data):
        """
        Set Earth position from JPL Horizons data.

        Args:
            position_data: Dictionary with Earth position components
        """
        self.xyz_earth = [
            FLOAT_TYPE(position_data['x']),
            FLOAT_TYPE(position_data['y']),
            FLOAT_TYPE(position_data['z'])
        ]

        # Store components for convenience
        self.xe = self.xyz_earth[0]
        self.ye = self.xyz_earth[1]
        self.ze = self.xyz_earth[2]

    def setup_coordinate_system(self, xarr, yarr, zarr, trueanarr, heliorbit, comet):
        """
        Set up coordinate system and observation parameters.

        Args:
            xarr: Comet x-positions array
            yarr: Comet y-positions array
            zarr: Comet z-positions array
            trueanarr: Comet true anomaly array
        """
        # Comet coordinates at observation date
        self.xcobs = FLOAT_TYPE(xarr[self.ntimes - 1])
        self.ycobs = FLOAT_TYPE(yarr[self.ntimes - 1])
        self.zcobs = FLOAT_TYPE(zarr[self.ntimes - 1])

        # Comet true anomaly at observation date
        self.thetacobs = FLOAT_TYPE(trueanarr[self.ntimes - 1])

        # Comet heliocentric distance at observation
        self.rcobs = FLOAT_TYPE(np.sqrt(self.xcobs**2 + self.ycobs**2 + self.zcobs**2))

        # Convert to equatorial coordinates
        self.rah = int(self.ra / 15.0)
        self.ram = FLOAT_TYPE((self.ra / 15.0 - self.rah) * 60.0)
        self.ras = FLOAT_TYPE((self.ram - int(self.ram)) * 60.0)
        self.ded = int(self.dec)
        self.dem = FLOAT_TYPE((self.dec - self.ded) * 60.0)
        self.des = FLOAT_TYPE((self.dem - int(self.dem)) * 60.0)

        # Calculate brightness scaling factors
        # Apply Schleicher phase function correction
        x = FLOAT_TYPE(self.phas_ang)
        # Use Horner's method for polynomial evaluation to reduce floating-point errors
        phase_corr = FLOAT_TYPE(10 ** (SCHLEICHER_COEFFS[0] +
                                       x * (SCHLEICHER_COEFFS[1] +
                                            x * (SCHLEICHER_COEFFS[2] +
                                                 x * (SCHLEICHER_COEFFS[3] +
                                                      x * (SCHLEICHER_COEFFS[4] +
                                                           x * (SCHLEICHER_COEFFS[5] +
                                                                x * SCHLEICHER_COEFFS[6])))))))

        self.pv = FLOAT_TYPE(self.pv0 * phase_corr)

        # Linear phase correction for nucleus
        self.pv_nuc = FLOAT_TYPE(self.pv0_nuc * 10 ** (-self.phase_coeff_nuc * self.phas_ang / 2.5))

        # Calculate brightness scaling constants
        self.cte_part = FLOAT_TYPE(self.pv * self.solid_angle /
                                  ((AUKM * 1.0e3) ** 2 * self.rcobs ** 2 * self.delta ** 2 * self.areapix))

        self.cte_nuc = FLOAT_TYPE(self.pv_nuc * self.solid_angle /
                                 ((AUKM * 1.0e3) ** 2 * self.rcobs ** 2 * self.delta ** 2 * self.areapix))

        # Calculate nucleus brightness
        self.brnucleus = FLOAT_TYPE(self.cte_nuc * self.rnucleus ** 2)

        # Convert to heliocentric orbit plane coordinates
        self.xoute, self.youte, self.zoute = he_to_hpo(self.xe, self.ye, self.ze, self.helio_matrix)

        # Cometocentric coordinates of Earth at observation date
        self.chitae = FLOAT_TYPE(self.xoute * np.cos(self.thetacobs) + self.youte * np.sin(self.thetacobs) - self.rcobs)
        self.etae = FLOAT_TYPE(self.xoute * np.sin(self.thetacobs) - self.youte * np.cos(self.thetacobs))
        self.gitae = FLOAT_TYPE(self.zoute)

        # Conversion parameters to sky plane (N-M)
        self.nmpar = FLOAT_TYPE(np.sqrt(self.etae ** 2 + self.gitae ** 2))
        self.nmpar1 = FLOAT_TYPE(self.nmpar / self.delta)
        self.nmpar2 = FLOAT_TYPE(self.chitae * self.etae / (self.nmpar * self.delta))
        self.nmpar3 = FLOAT_TYPE(self.chitae * self.gitae / (self.nmpar * self.delta))
        self.nmpar4 = FLOAT_TYPE(self.gitae / self.nmpar)
        self.nmpar5 = FLOAT_TYPE(self.etae / self.nmpar)

        self.nmpar6 = FLOAT_TYPE(self.chitae / self.delta)
        self.nmpar7 = FLOAT_TYPE(self.etae / self.delta)
        self.nmpar8 = FLOAT_TYPE(self.gitae / self.delta)

        # Pass parameters to the heliorbit object
        heliorbit.set_matrices(self.helio_matrix, None)
        heliorbit.set_params_nm(
            self.delta, self.nmpar, self.nmpar1, self.nmpar2, self.nmpar3,
            self.nmpar4, self.nmpar5, self.nmpar6, self.nmpar7, self.nmpar8)

        # Observation time since perihelion passage in days
        self.tc = FLOAT_TYPE(self.end_jd - self.per_jd)

        # Size step
        self.r1log = FLOAT_TYPE(np.log10(self.radiomin[0]))
        self.r2log = FLOAT_TYPE(np.log10(self.radiomax[0]))
        self.steprlog = FLOAT_TYPE((self.r2log - self.r1log) / float(self.nsizes))

        # Add missing output information
        print(f"  Helioc. coor. of the Earth at obs. date:   {self.xe:8.7f}  {self.ye:8.7f}  {self.ze:8.7f}")
        print(f"  TARGET RA (h,m,s):    {self.rah:2d} {int(self.ram):2d} {self.ras:5.2f}")
        print(f"  TARGET DEC(d,m,s):   {self.ded:2d} {int(self.dem):2d} {self.des:5.2f}")
        print(
            f"  Phase angle = {self.phas_ang:7.3f} Pos. angle = {self.psang:7.3f} Plane angle= {self.plang:7.3f} (deg)")
        print(f"  Comet elements: EC= {comet.ec:9.6f} QR= {comet.qr:9.6f} au    TP= {comet.per_jd:12.5f}")
        print(f"  Comet elements: NODE= {comet.om:7.3f} W= {comet.wper:7.3f} INC= {comet.inc:7.3f} (deg)")
        print(f"  Target heliocentric distance= {self.rcobs:10.3f} au  Distance to observer= {self.delta:10.3f} au")