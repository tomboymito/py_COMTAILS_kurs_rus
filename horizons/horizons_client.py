"""
JPL Horizons API client for COMTAILS simulation.

This module provides a client for downloading ephemeris data from
the JPL Horizons system for comets and observation points.
"""
import os
import re
import requests

from constants import FLOAT_TYPE


class HorizonsClient:
    """
    Client for interacting with the JPL Horizons API.

    This class handles the communication with the JPL Horizons system
    to fetch orbital elements and position data for solar system bodies.
    """

    def __init__(self, base_url="https://ssd.jpl.nasa.gov/api/horizons_file.api"):
        """
        Initialize the Horizons client.

        Args:
            base_url: URL for the JPL Horizons API
        """
        self.base_url = base_url
        self.output_dir = "output"

        # Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

    def get_earth_position(self, jd):
        """
        Get Earth position at the specified Julian date.

        Args:
            jd: Julian date

        Returns:
            dict: Dictionary with Earth position components
        """
        # Create ephemeris request
        query_data = self._create_observer_request(jd)

        # Write request to file for reference
        with open(f'{self.output_dir}/obs_point.txt', 'w') as f:
            f.write(query_data)

        # Send request to JPL
        response = self._send_request(query_data)

        # Parse Earth position
        earth_position = self._parse_earth_position(response)

        return earth_position

    def get_comet_data(self, comet_id, jd):
        """
        Get comet data at the specified Julian date.

        Args:
            comet_id: JPL Horizons record ID for the comet
            jd: Julian date

        Returns:
            dict: Dictionary with comet orbital elements and other data
        """
        # Create ephemeris request
        query_data = self._create_comet_request(comet_id, jd)

        # Write request to file for reference
        with open(f'{self.output_dir}/comet_ra_dec_elements.txt', 'w') as f:
            f.write(query_data)

        # Send request to JPL
        response = self._send_request(query_data)

        # Parse comet data
        comet_data = self._parse_comet_data(response)

        return comet_data

    def _create_observer_request(self, jd):
        """
        Create JPL Horizons request for Earth position.

        Args:
            jd: Julian date

        Returns:
            str: Request string
        """
        request = "!$$SOF\n"
        request += "COMMAND='399'\n"  # Set to Earth
        request += "OBJ_DATA='NO'\n"
        request += "MAKE_EPHEM='YES'\n"
        request += "REF_PLANE='ECLIPTIC'\n"
        request += "TABLE_TYPE='VECTORS'\n"
        request += "CENTER='500@10'\n"  # Center=Sun
        request += "OUT_UNITS='AU-D'\n"
        request += "VEC_TABLE='2'\n"
        request += f"TLIST='{jd:.6f}'\n"

        return request

    def _create_comet_request(self, comet_id, jd):
        """
        Create JPL Horizons request for comet data.

        Args:
            comet_id: JPL Horizons record ID for the comet
            jd: Julian date

        Returns:
            str: Request string
        """
        request = "!$$SOF\n"
        request += f"COMMAND='{comet_id}';\n"
        request += "OBJ_DATA='NO'\n"
        request += "EPHEM_TYPE='OBS'\n"
        request += "REF_PLANE='ECLIPTIC'\n"
        request += "ANG_FORMAT='DEG'\n"
        request += "CENTER='399'\n"  # geocentric
        request += "OUT_UNITS='AU-D'\n"
        request += "CAL_FORMAT='JD'\n"
        request += "QUANTITIES='1,20,27,28,43'\n"  # RA,DEC,PsAng,PlAng,Phase
        request += f"TLIST='{jd:.6f}'\n"

        return request

    def _send_request(self, query_data):
        """
        Send request to JPL Horizons API.

        Args:
            query_data: Request string

        Returns:
            str: Response text
        """
        params = {"format": "text", "input": query_data}

        try:
            response = requests.post(self.base_url, data=params)

            if response.status_code != 200:
                print(f"Error getting data from JPL: {response.status_code}")
                return ""

            return response.text

        except Exception as e:
            print(f"Error in JPL request: {e}")
            return ""

    def _parse_earth_position(self, response):
        """
        Parse Earth position from JPL Horizons response.

        Args:
            response: Response text from JPL

        Returns:
            dict: Dictionary with Earth position components
        """
        # Default values
        earth_position = {
            'x': FLOAT_TYPE(0.0),
            'y': FLOAT_TYPE(0.0),
            'z': FLOAT_TYPE(0.0)
        }

        # Parse response
        lines = response.splitlines()
        earth_position_found = False

        for i, line in enumerate(lines):
            # Look for the vector data section
            if "$$SOE" in line:
                for j in range(i + 1, min(i + 10, len(lines))):
                    vector_line = lines[j]

                    # Check if this is an X,Y,Z position line
                    if "X =" in vector_line and "Y =" in vector_line and "Z =" in vector_line:
                        # Extract values with labels X =, Y =, Z =
                        try:
                            x_part = vector_line.split("X =")[1].split("Y =")[0].strip()
                            y_part = vector_line.split("Y =")[1].split("Z =")[0].strip()
                            z_part = vector_line.split("Z =")[1].strip()

                            earth_position['x'] = FLOAT_TYPE(float(x_part))
                            earth_position['y'] = FLOAT_TYPE(float(y_part))
                            earth_position['z'] = FLOAT_TYPE(float(z_part))

                            print(f"Earth position found: {earth_position}")
                            earth_position_found = True
                            break
                        except (ValueError, IndexError) as e:
                            print(f"Failed to parse Earth position: {e}")
                            continue
                break

        if not earth_position_found:
            print("WARNING: Could not find valid Earth position in JPL response")

        return earth_position

    def _parse_comet_data(self, response):
        """
        Parse comet data from JPL Horizons response.

        Args:
            response: Response text from JPL

        Returns:
            dict: Dictionary with comet orbital elements and other data
        """
        # Default values for orbital elements
        comet_data = {
            'ec': FLOAT_TYPE(0.5),  # eccentricity
            'qr': FLOAT_TYPE(1.0),  # perihelion distance (AU)
            'per_jd': FLOAT_TYPE(0.0),  # perihelion time
            'om': FLOAT_TYPE(0.0),  # longitude of ascending node (deg)
            'wper': FLOAT_TYPE(0.0),  # argument of perihelion (deg)
            'inc': FLOAT_TYPE(0.0),  # inclination (deg)
            'ra': FLOAT_TYPE(0.0),  # right ascension (deg)
            'dec': FLOAT_TYPE(0.0),  # declination (deg)
            'delta': FLOAT_TYPE(1.0),  # distance to observer (AU)
            'deldot': FLOAT_TYPE(0.0),  # rate of change of distance
            'psang': FLOAT_TYPE(0.0),  # position angle
            'psamv': FLOAT_TYPE(0.0),  # position angle of heliocentric velocity
            'plang': FLOAT_TYPE(0.0),  # phase angle of orbital plane
            'phas_ang': FLOAT_TYPE(0.0)  # phase angle
        }

        # Parse response
        lines = response.splitlines()

        for i, line in enumerate(lines):
            if "Target body name" in line:
                print(" " + line)

            # Look for the EPOCH section with orbital elements
            if "EC=" in line:
                # This line likely contains EC, QR, and TP
                ec_match = re.search(r'EC= ([0-9.]+)', line)
                if ec_match:
                    comet_data['ec'] = FLOAT_TYPE(float(ec_match.group(1)))

                qr_match = re.search(r'QR= ([0-9.]+)', line)
                if qr_match:
                    comet_data['qr'] = FLOAT_TYPE(float(qr_match.group(1)))

                tp_match = re.search(r'TP= ([0-9.]+)', line)
                if tp_match:
                    comet_data['per_jd'] = FLOAT_TYPE(float(tp_match.group(1)))

            if "OM=" in line:
                # This line likely contains OM, W, and IN
                om_match = re.search(r'OM= ([0-9.]+)', line)
                if om_match:
                    comet_data['om'] = FLOAT_TYPE(float(om_match.group(1)))

                # Try different patterns for argument of perihelion
                w_match = re.search(r'W = ([0-9.]+)', line)
                if w_match:
                    comet_data['wper'] = FLOAT_TYPE(float(w_match.group(1)))
                else:
                    # Alternative pattern
                    w_match = re.search(r'W= ([0-9.]+)', line)
                    if w_match:
                        comet_data['wper'] = FLOAT_TYPE(float(w_match.group(1)))

                inc_match = re.search(r'IN= ([0-9.]+)', line)
                if inc_match:
                    comet_data['inc'] = FLOAT_TYPE(float(inc_match.group(1)))

            # Extract astrometric data
            if "$SOE" in line and i + 1 < len(lines):
                # Get astrometric data (next line after SOE)
                next_line = lines[i + 1]

                try:
                    parts = next_line[22:].split()
                    if len(parts) >= 8:
                        comet_data['ra'] = FLOAT_TYPE(float(parts[0]))
                        comet_data['dec'] = FLOAT_TYPE(float(parts[1]))
                        comet_data['delta'] = FLOAT_TYPE(float(parts[2]))
                        comet_data['deldot'] = FLOAT_TYPE(float(parts[3]))
                        comet_data['psang'] = FLOAT_TYPE(float(parts[4]))
                        comet_data['psamv'] = FLOAT_TYPE(float(parts[5]))
                        comet_data['plang'] = FLOAT_TYPE(float(parts[6]))
                        comet_data['phas_ang'] = FLOAT_TYPE(float(parts[7]))
                except Exception as e:
                    print(f"Failed to parse astrometric data: {e}")
                    # Default values are already set

                break

        return comet_data
