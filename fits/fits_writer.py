"""
FITS file writing module for COMTAILS simulation.

This module provides utilities for writing FITS images from simulation results.
"""
import os
import numpy as np
import astropy.io.fits as fits

class FitsWriter:
    """
    Class for writing FITS files from simulation data.

    This class provides methods to write simulation results to FITS files
    with appropriate header information.
    """

    def __init__(self):
        """Initialize the FITS writer."""
        # Ensure output directory exists
        os.makedirs('output', exist_ok=True)

    def write_fits_image(self, outimage, imagefit, object_name, start_jd, end_jd,
                         grdsiz, ntotmc, swap_axis_and_subtract_1_1=True):
        """
        Write a FITS image to disk.

        Args:
            outimage: Path to output FITS file
            imagefit: Image array to write
            object_name: Name of the object (comet)
            start_jd: Start Julian date
            end_jd: End Julian date (observation)
            grdsiz: Grid size in km/pixel
            ntotmc: Total number of Monte Carlo events
            swap_axis_and_subtract_1_1: Whether to swap x and y axes and shift entire matrix by (-1,-1)
                                      (default: True to match FORTRAN results)
        """
        # Make a copy to avoid modifying the original data
        data_to_write = imagefit.copy()

        # Swap axes if requested
        if swap_axis_and_subtract_1_1:
            data_to_write = data_to_write.transpose()
            # Create a new array with shift applied - pad with zeros
            # This effectively shifts the entire image by (-1,-1)
            shifted_data = np.zeros_like(data_to_write)
            shifted_data[:-1, :-1] = data_to_write[1:, 1:]
            data_to_write = shifted_data

        # Create a new FITS primary array
        hdu = fits.PrimaryHDU(data_to_write)

        # Add header information
        hdu.header['OBJECT'] = object_name
        hdu.header['START_TM'] = start_jd
        hdu.header['OBS_TIME'] = end_jd
        hdu.header['PIXSIZ'] = grdsiz
        hdu.header['MC-Event'] = ntotmc
        hdu.header['SWAPAXES'] = swap_axis_and_subtract_1_1

        # Write to disk
        hdu.writeto(outimage, overwrite=True)
        print(f"Wrote FITS file: {outimage}")
