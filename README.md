# COMTAILS - COMetary dust TAIL Simulator

COMTAILS is a comprehensive Python implementation of a comet dust tail simulation program. It improves upon the original FORTRAN 77 code developed by Fernando Moreno (IAA-CSIC) with a cleaner object-oriented design, enhanced numerical stability, and improved visualization capabilities. The porting of the FORTRAN serial version to this parallel Python version was performed by Rafael Morales and Nicolás Robles (IAA-CSIC).
The original reposoritory with FORTRAN versions is: https://github.com/FernandoMorenoDanvila/COMTAILS/tree/FORTRAN_SERIAL (FORTRAN serial)  and https://github.com/FernandoMorenoDanvila/COMTAILS/tree/FORTRAN_PARALLEL (FORTRAN MPI parallel version).

## Description

COMTAILS generates realistic simulations of comet dust tails by modeling the dynamics of dust particles under the influence of solar radiation pressure and solar gravity. It produces high-quality images that can be compared with actual observations of comets.

### Key Features

- Monte Carlo simulation of dust particle dynamics 
- High-precision orbital calculations using Kepler's equation solvers
- Support for different ejection models (isotropic, sunward, active areas)
- Integration with JPL Horizons for accurate ephemerides
- Star field overlay from Gaia EDR3
- FITS image output for scientific analysis
- Particle visualization with pygame
- Multiprocessing support for improved performance
- Calculation of Afrho and magnitude parameters

## Installation

### Prerequisites

- Python 3.8 or later
- NumPy
- AstroPy
- PyGame (for visualization)
- Requests (for JPL Horizons API access)

### Installation Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/username/comtails.git
   cd comtails
   ```

2. Install dependencies:
   ```bash
   pip install numpy astropy pygame requests
   ```

## Usage

### Basic Usage

Run a simulation with default parameters:

```bash
python main.py
```

### Custom Configuration

Specify custom configuration files:

```bash
python main.py --input-dir custom_inputs --config my_config.dat --dust-profile my_profile.dat
```

### Configuration Files

COMTAILS requires two main configuration files:

1. **Main configuration file** (default: `TAIL_INPUTS.dat`):
   - Contains comet parameters, observation setup, and simulation controls
   - Specifies date range, image grid size, and dust physical properties

2. **Dust loss rate profile** (default: `dmdt_vel_power_rmin_rmax.dat`):
   - Defines dust production rates, velocity factors, size distribution power-law indices, and particle size ranges vs. time

## Output Files

The simulation produces the following output files in the `output` directory:

- `tail_sdu.fits`: Dust tail brightness image (solar disk intensity units)
- `tail_mag.fits`: Dust tail brightness image (mag/arcsec²)
- `OPT_DEPTH.fits`: Optical depth map
- `afrho.dat`: Afrho parameter vs. time
- `dust_particles.png`: Visualization of particle positions (if enabled)
- `dustlossrate.dat`: Dust loss rate vs. time
- Additional files based on configuration settings

## Code Structure

- `main.py`: Entry point for running the simulation
- `simulation.py`: Main simulation controller
- `config.py`: Configuration management
- `constants.py`: Physical and mathematical constants
- `comet.py`: Comet model and orbital properties
- `dust_tail.py`: Dust tail simulation with Monte Carlo methods
- `heliorbit.py`: Heliocentric orbit calculations
- `orbit_solver.py`: Kepler's equation solvers
- `star_field.py`: Star field generation
- `plot_handler.py`: Visualization utilities
- `horizons_client.py`: JPL Horizons API client
- `fits_writer.py`: FITS image output

## License

See the LICENSE file for details (MIT license).

## Citation

If you use this code in your research, please write in  your paper that:

``The model results are based on a python implementation, performed by Rafael Morales and Nicolás Robles  of the Instituto de Astrofísica de Andalucía, from the original FORTRAN serial code written by Fernando Moreno (Moreno, 2025)``

and include in your reference list:

``Moreno, F. (2025). COMetary dust TAIL Simulator (COMTAILS): A computer code to generate comet dust tail brightness images
Astronomy and Astrophysics, Volume:695(2025), Article:A263.``


## Acknowledgments

If you use this code in your research, please write in the acknowledgements of your papert:

``The model results are based on a python implementation, performed by Rafael Morales and Nicolás Robles of
the Instituto de Astrofísica de Andalucía, from the original FORTRAN serial code written by Fernando Moreno (Moreno, 2025)``

and give also acknowledgements to:

- NASA/JPL-Caltech for the Horizons ephemeris system
- ESA/Gaia for the star catalog data
