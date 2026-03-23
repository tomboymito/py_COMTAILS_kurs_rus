"""
Main simulation controller for COMTAILS.

This module provides the core simulation controller that orchestrates the
comet dust tail simulation. It improves upon the original COMTAILS.for Fortran 77 code
by Fernando Moreno IAA-CSIC by implementing a cleaner object-oriented design.
"""
import time
import numpy as np

from config import SimulationConfig
from fits.fits_writer import FitsWriter
from horizons.horizons_client import HorizonsClient
from models.comet import Comet
from models.dust_tail import DustTail
from visualization.star_field import StarField
from utils.coordinate_transforms import TransformationFactory
from constants import FLOAT_TYPE


class SimulationController:
    """
    Main controller for the COMTAILS simulation.

    This class coordinates the entire simulation pipeline from reading inputs to
    producing the final outputs.
    """

    def __init__(self):
        """Initialize the simulation controller."""
        self.config = None
        self.comet = None
        self.dust_tail = None
        self.horizons_client = None
        self.fits_writer = None
        self.star_field = None
        self.plot_handler = None
        self.transform_factory = None

        # Performance tracking
        self.start_time = None
        self.end_time = None

    def run(self, config_files):
        """
        Run the full simulation pipeline.

        Args:
            config_files: List of configuration file paths

        Returns:
            dict: Simulation results summary
        """
        # Start CPU timer
        self.start_time = time.time()

        # Initialize components
        self._init_components(config_files)

        # Process simulation steps
        self._setup_time_array()
        self._fetch_comet_coordinates()
        self._setup_image_grid()
        self._setup_coordinate_system()
        self._process_star_field()
        self._build_dust_tail()
        self._finalize_image()

        # Write results
        results = self._write_results()

        # End CPU timer
        self.end_time = time.time()
        elapsed_minutes = (self.end_time - self.start_time) / 60.0
        print(f" Elapsed CPU time: {elapsed_minutes:12.5f} minutes")

        return results

    def _init_components(self, config_files):
        """Initialize all simulation components."""
        # Read configuration
        self.config = SimulationConfig()
        self.config.read_inputs(config_files)

        # Initialize models
        self.comet = Comet(self.config)
        self.dust_tail = DustTail(self.config)

        # Initialize I/O components
        self.horizons_client = HorizonsClient()
        self.fits_writer = FitsWriter()

        # Initialize visualization components
        if self.config.istar == 1:
            self.star_field = StarField(self.config)

        # Initialize utilities
        self.transform_factory = TransformationFactory()

    def _setup_time_array(self):
        """Set up time array for the simulation."""
        self.config.setup_time_array()

    def _fetch_comet_coordinates(self):
        """Fetch comet coordinates and orbital elements from JPL Horizons."""
        print("Downloading ephemeris from JPL-Horizons...")

        # Get Earth position
        earth_data = self.horizons_client.get_earth_position(self.config.end_jd)
        self.config.set_earth_position(earth_data)

        # Get comet coordinates and elements
        comet_data = self.horizons_client.get_comet_data(
            self.config.recid,
            self.config.end_jd
        )
        self.comet.set_orbital_elements(comet_data)

        # Calculate heliocentric matrix
        self.config.helio_matrix = self.transform_factory.set_helio_matrix(
            self.comet.inc_rad,
            self.comet.om_rad,
            self.comet.wper_rad
        )

        # Compute positions and velocities for all times
        self.comet.compute_positions_and_velocities(
            self.config.times,
            self.config.helio_matrix
        )

        # Print orbit information
        self.comet.print_orbit_info(
            self.config.start_jd,
            self.config.end_jd,
            self.config.per_jd
        )

    def _setup_image_grid(self):
        """Set up the image grid for the dust tail simulation."""
        self.plot_handler = self.config.setup_image_grid(self.comet)

    def _setup_coordinate_system(self):
        """Set up coordinate systems and observation parameters."""
        self.config.setup_coordinate_system(
            self.comet.xarr,
            self.comet.yarr,
            self.comet.zarr,
            self.comet.trueanarr,
            self.dust_tail.heliorbit,
            self.comet
        )

    def _process_star_field(self):
        """Process star field within the image."""
        if self.config.istar == 1 and self.star_field:
            self.star_field.download_star_field()
            self.star_field.process_star_field()
            self.config.flux_star = self.star_field.get_flux_array()

    def _build_dust_tail(self):
        """Build the dust tail through Monte Carlo simulation."""
        # Fixed formatting to match expected output
        print("Building up dust tail ...")
        self.dust_tail.build(
            self.comet,
            self.config,
            self.plot_handler
        )

        # Add the missing output line for total dust mass
        print(f"  Total dust mass ejected= {self.dust_tail.total_mass:.3E} kg")

    def _finalize_image(self):
        """Finalize the dust tail image and apply effects."""
        # Add nucleus contribution, attenuated by optical thickness
        self.config.flux[self.config.inuc, self.config.jnuc] += FLOAT_TYPE(
            self.config.brnucleus * np.exp(-self.dust_tail.opt_depth_nuc)
        )

        # Attenuate star flux by optical depth of dust tail
        for ii in range(self.config.nx):
            for jj in range(self.config.ny):
                self.config.flux_star[ii, jj] = FLOAT_TYPE(
                    self.config.flux_star[ii, jj] *
                    np.exp(-self.dust_tail.opt_depth[ii, jj])
                )

        # Add star field to total flux
        star_flux_sum = np.sum(self.config.flux_star)
        print(
            f"Star field contains {np.count_nonzero(self.config.flux_star)} non-zero pixels, max flux: {np.max(self.config.flux_star)}")
        print(f"Adding total star flux: {star_flux_sum}")
        self.config.flux += self.config.flux_star

        # Apply convolution if requested
        if self.config.iconv == 1:
            self.dust_tail.apply_convolution(self.config.flux, self.config.sfwhm)

    def _write_results(self):
        """Write simulation results to files."""
        # Compute Afrho and magnitude
        afrho, afrho_0, mag = self.dust_tail.calculate_afrho_mag(
            self.config.flux,
            self.config.inuc,
            self.config.jnuc,
            self.config.grdsiz,
            self.config.cte_mag,
            self.config.magsun,
            self.config.scale,
            self.config.delta,
            self.config.rcobs,
            self.config.rho_ap
        )

        # Fixed formatting to match expected output
        print(f"  Aperture(km)= {self.config.rho_ap:10.2f} Afrho(m)= {afrho:8.3f} Mag= {mag:8.3f}")

        # Write output files
        self.fits_writer.write_fits_image(
            'output/tail_sdu.fits',
            self.config.flux,
            self.config.recid,
            self.config.start_jd,
            self.config.end_jd,
            self.config.grdsiz,
            self.config.ntotmc
        )

        # Write magnitude FITS file
        flux_mag = self._calculate_magnitude_image()
        self.fits_writer.write_fits_image(
            'output/tail_mag.fits',
            flux_mag,
            self.config.recid,
            self.config.start_jd,
            self.config.end_jd,
            self.config.grdsiz,
            self.config.ntotmc
        )

        # Write optical depth FITS file
        self.fits_writer.write_fits_image(
            'output/OPT_DEPTH.fits',
            self.dust_tail.opt_depth,
            self.config.recid,
            self.config.start_jd,
            self.config.end_jd,
            self.config.grdsiz,
            self.config.ntotmc
        )

        # Save particle plot if enabled
        if self.config.igrapho == 1 and self.plot_handler and self.plot_handler.available:
            self.plot_handler.save_image()
            self.plot_handler.close()
            print("Particle plot generated successfully")

        # Write Afrho to file
        with open("output/afrho.dat", "a") as f:
            f.write(f"{self.config.end_jd - self.config.per_jd:20.12e} "
                    f"{afrho:15.8e} {afrho_0:15.8e} {mag:15.8e}\n")

        # Return results summary
        return {
            "afrho": afrho,
            "afrho_0": afrho_0,
            "mag": mag,
            "total_dust_mass": self.dust_tail.total_mass
        }

    def _calculate_magnitude_image(self):
        """Calculate image in mag/arcsec^2 units."""
        flux_mag = np.zeros((self.config.nx, self.config.ny), dtype=FLOAT_TYPE)

        for ii in range(self.config.nx):
            for jj in range(self.config.ny):
                if self.config.flux[ii, jj] > 0.0:
                    flux_mag[ii, jj] = FLOAT_TYPE(
                        self.config.cte_mag + self.config.magsun -
                        2.5 * np.log10(self.config.flux[ii, jj])
                    )
                else:
                    flux_mag[ii, jj] = FLOAT_TYPE(100.0)  # mag/arcsec^2 (i.e., ~zero flux)

        return flux_mag

    @staticmethod
    def validate_results(results_file="output/afrho.dat", expected_afrho=10.5,
                         expected_mag=8.07, tolerance=0.1):
        """
        Validate the results of a COMTAILS simulation.

        Args:
            results_file: Path to results file
            expected_afrho: Expected Afrho value
            expected_mag: Expected magnitude
            tolerance: Acceptable relative difference

        Returns:
            bool: True if results are within tolerance, False otherwise
        """
        try:
            # FIX: Pass the file path as a string, not an object
            with open(results_file, 'r') as f:
                last_line = f.readlines()[-1]
                parts = last_line.split()
                afrho = float(parts[1])
                mag = float(parts[3])

            afrho_diff = abs(afrho - expected_afrho) / expected_afrho
            mag_diff = abs(mag - expected_mag) / expected_mag

            if afrho_diff > tolerance or mag_diff > tolerance:
                print(f"Validation failed: Afrho={afrho} (expected {expected_afrho}), "
                      f"Mag={mag} (expected {expected_mag})")
                return False
            else:
                print(f"Validation passed: Afrho={afrho}, Mag={mag}")
                return True
        except Exception as e:
            print(f"Validation error: {e}")
            return False
