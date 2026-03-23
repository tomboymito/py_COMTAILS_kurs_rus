"""
Star field handling module for COMTAILS simulation.

This module provides functionality for downloading star field data
from Gaia EDR3 and processing it for inclusion in comet images.
"""
import os
import requests
import numpy as np

from constants import FLOAT_TYPE
from utils.coordinate_transforms import std_coor


class StarField:
    """
    Class for downloading and processing star field data.

    This class handles fetching star positions from catalogs
    and converting them to the image coordinate system.
    """

    def __init__(self, config):
        """
        Initialize the star field handler.

        Args:
            config: SimulationConfig object
        """
        self.config = config
        self.output_dir = "output"
        self.stars = []
        self.flux_array = np.zeros((config.nx, config.ny), dtype=FLOAT_TYPE)
        self.star_count = 0

        # Create output directory
        os.makedirs(self.output_dir, exist_ok=True)

    def download_star_field(self, apply_filtering=False, max_stars=10000):
        """
        Download star field data from Gaia EDR3.

        Args:
            apply_filtering: Whether to filter stars by magnitude and count
            max_stars: Maximum number of stars to return

        Returns:
            bool: True if successful, False otherwise
        """
        # Format coordinates for URL - include RA and DEC in format expected by API
        ra_str = f"{self.config.ra:.3f}"
        de_str = f"{self.config.dec:.3f}"

        # Include sign for declination formatting
        sign = '+' if self.config.dec >= 0 else ''

        # Construct URL for Gaia EDR3 data
        base_url = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"

        # Use the proper format for the objstr parameter (RA DEC)
        # The API expects decimal degrees without the 'd' suffix
        params = {
            "catalog": "gaia_edr3_source",
            "spatial": "cone",
            "objstr": f"{ra_str} {sign}{de_str}",
            "radius": "3600",  # 1 degree in arcseconds
            "outfmt": "1",
            "selcols": "ra,dec,phot_g_mean_mag,phot_rp_mean_mag,phot_bp_mean_mag"
        }

        # Debug: Print the request parameters
        print(f"Query coordinates: RA={ra_str}, DEC={sign}{de_str}")
        print(f"Magnitude limit: {self.config.maglim}")

        # Download data
        try:
            print("Downloading star field coordinate/mag table...")
            response = requests.get(base_url, params=params)

            if response.status_code != 200:
                print(f"Error downloading star data: HTTP {response.status_code}")
                return False

            # Check for server error in response
            if 'stat="ERROR"' in response.text:
                print(f"Server returned an error in response")
                print(response.text[:500])  # Print first part of response for debugging
                return False

            # Save raw data to file
            star_data_file = f'{self.output_dir}/star.dat'
            with open(star_data_file, 'w') as f:
                f.write(response.text)

            print(f"Star data saved to {star_data_file}")

            # Check if the file contains actual data
            with open(star_data_file, 'r') as f:
                content = f.read()
                if '|' not in content:
                    print("Warning: Downloaded file does not contain expected column separator '|'")
                    print("First 100 characters of response:")
                    print(content[:100])
                else:
                    # Count data lines (very rough estimate)
                    data_lines = sum(1 for line in content.splitlines() if '|' in line)
                    print(f"Downloaded approximately {data_lines} data lines")

            return True

        except Exception as e:
            print(f"Error in download_star_field: {e}")
            return False

    def process_star_field(self, apply_filtering=False):
        """
        Process downloaded star field data.

        Args:
            apply_filtering: Whether to filter stars by magnitude

        Returns:
            int: Number of stars processed
        """
        star_data_file = f'{self.output_dir}/star.dat'
        output_file = f'{self.output_dir}/starpos.dat'

        # Reset flux array
        self.flux_array = np.zeros((self.config.nx, self.config.ny), dtype=FLOAT_TYPE)

        try:
            # Open output file
            with open(output_file, 'w') as star_file:
                print(f"Processing star data from {star_data_file}")

                # Read the entire file to analyze structure
                with open(star_data_file, 'r') as f:
                    lines = f.readlines()

                # Find the header-data separator line (contains column names)
                header_end = 0
                for i, line in enumerate(lines):
                    if '|' in line and 'ra' in line.lower() and 'dec' in line.lower():
                        header_end = i
                        break

                # Skip 4 more lines to get to actual data (column names, types, units, null line)
                data_start = header_end + 4

                # Process each data line
                self.star_count = 0
                processed_count = 0  # Count all processed stars regardless of magnitude

                for i in range(data_start, len(lines)):
                    line = lines[i].strip()
                    if not line or line.startswith('\\'):
                        continue

                    # Split the line by whitespace (multiple spaces treated as one)
                    parts = [p for p in line.split() if p]

                    if len(parts) < 7:  # Need at least RA, DEC, and magnitudes
                        continue

                    try:
                        # Extract star data from columns
                        ra_star = float(parts[0])
                        de_star = float(parts[1])

                        # Handle 'null' values in magnitudes
                        if parts[4].lower() == 'null' or parts[5].lower() == 'null' or parts[6].lower() == 'null':
                            continue

                        # The magnitudes are in specific columns
                        gmag = float(parts[4])  # G magnitude
                        rpmag = float(parts[5])  # RP magnitude
                        bpmag = float(parts[6])  # BP magnitude

                        processed_count += 1

                        # Skip stars outside magnitude limit early
                        if apply_filtering and gmag > self.config.maglim:
                            continue

                        # Transform to standard coordinates
                        xtemp, ytemp = std_coor(
                            self.config.ra,
                            self.config.dec,
                            ra_star,
                            de_star
                        )

                        # Convert to pixel coordinates
                        posx = int(xtemp * self.config.arcsec_rad / self.config.scale) + self.config.inuc
                        posy = int(ytemp * self.config.arcsec_rad / self.config.scale) + self.config.jnuc

                        # Invert coordinates to match FORTRAN convention
                        posx = self.config.nx - posx
                        posy = self.config.ny - posy

                        # Convert to R-Cousins magnitude using Gaia colors
                        # Formula from Jordi et al. 2010, A&A, 523, A48
                        xtemp = bpmag - rpmag
                        ytemp = (0.02275 + 0.3691 * xtemp - 0.1243 * xtemp ** 2 -
                                0.01396 * xtemp ** 3 + 0.003775 * xtemp ** 4)
                        rmag = gmag - ytemp

                        starmag = rmag
                        starflux = 10 ** (-0.4 * (starmag + 10.925))

                        # Check if star is within image boundaries
                        if (posx >= 0 and posx < self.config.nx and
                                posy >= 0 and posy < self.config.ny):

                            # Apply magnitude filter
                            if starmag < self.config.maglim:
                                # Write to output file
                                star_file.write(
                                    f"{posx:5d} {posy:5d} {starmag:12.4e} {starflux:12.4e}\n")

                                # Add flux to array
                                self.flux_array[posx, posy] += starflux

                                # Store star information
                                self.stars.append({
                                    'ra': ra_star,
                                    'dec': de_star,
                                    'x': posx,
                                    'y': posy,
                                    'mag': starmag,
                                    'flux': starflux
                                })

                                self.star_count += 1

                    except Exception as e:
                        # Just skip problematic lines silently
                        continue

                print(f"Processed {processed_count} stars, {self.star_count} added within magnitude limit {self.config.maglim}")
                return self.star_count

        except Exception as e:
            print(f"Error processing star field: {e}")
            return 0

    def get_flux_array(self):
        """
        Get the star flux array.

        Returns:
            ndarray: Star flux array
        """
        return self.flux_array
