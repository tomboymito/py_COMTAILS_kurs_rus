import numpy as np
import multiprocessing
import os

from constants import (
    FLOAT_TYPE, PI, HALFPI, TWOPI, FOURPI, TORAD,
    CTEVEL, QP, SCHLEICHER_COEFFS, AUKM
)
from orbital.heliorbit import HelioOrbit
from utils.numerical import powerint

class DustTail:
    """
    Class for modeling a comet dust tail.

    This class handles the dust tail simulation through Monte Carlo methods.
    """

    def __init__(self, config):
        """
        Initialize dust tail object.

        Args:
            config: SimulationConfig object containing configuration parameters
        """
        # Initialize arrays
        self.opt_depth = np.zeros((config.nx, config.ny), dtype=FLOAT_TYPE)
        self.opt_depth_nuc = FLOAT_TYPE(0.0)

        # Initialize orbit handler
        self.heliorbit = HelioOrbit()

        # Initialize total mass
        self.total_mass = FLOAT_TYPE(0.0)

    def build(self, comet, config, plot_handler=None):
        """
        Build the dust tail through Monte Carlo simulation.

        Args:
            comet: Comet object
            config: SimulationConfig object
            plot_handler: Optional PlotHandler for visualization
        """
        # Reset total mass at the beginning
        self.total_mass = FLOAT_TYPE(0.0)

        # Determine number of processes
        n_proc = multiprocessing.cpu_count()

        if n_proc > 1 and config.ntimes > 1:
            # Parallel version
            chunk_size = config.ntimes // n_proc
            ranges = [(i * chunk_size, (i + 1) * chunk_size if i < n_proc - 1 else config.ntimes)
                      for i in range(n_proc)]

            with multiprocessing.Pool(n_proc) as pool:
                results = pool.starmap(
                    dust_tail_worker,
                    [(start, end, comet, config) for start, end in ranges]
                )

            # Combine results
            config.flux = np.sum([res['flux_local'] for res in results], axis=0)
            self.opt_depth = np.sum([res['opt_depth_local'] for res in results], axis=0)
            self.total_mass = sum(res['total_mass_local'] for res in results)
            self.opt_depth_nuc = sum(res['opt_depth_nuc_local'] for res in results)

            # Collect and sort dustloss data
            dustloss_data = [item for res in results for item in res['dustloss_data']]
            dustloss_data.sort(key=lambda x: x[0])
            with open('output/dustlossrate.dat', 'w') as f:
                for time_rel, rc_ejec, dmdt in dustloss_data:
                    f.write(f"    {time_rel:10.5E}     {rc_ejec:10.5E}     {dmdt:10.5E}\n")

            # Handle particle positions output
            if config.iprn == 1:
                # Each process wrote to its own file; optionally combine them
                pass  # Leave as separate files for simplicity

            # Handle plotting
            if config.igrapho == 1 and plot_handler is not None and plot_handler.available:
                for res in results:
                    for npar, mpar in res['particles_to_plot']:
                        plot_handler.add_particle(npar, mpar)

            # Handle subsolar data if applicable
            if config.iejec_mode == 3:
                subsolar_data = [item for res in results for item in res['subsolar_data']]
                subsolar_data.sort(key=lambda x: x[0])
                with open('output/subsolar.dat', 'w') as f:
                    for time_rel, rc_ejec, theta, subsolat in subsolar_data:
                        f.write(f"{time_rel:15.5e} {rc_ejec:15.5e} {theta:15.5e} {subsolat:15.5e}\n")

        else:
            # Sequential version (original code)
            for itime in range(config.ntimes - 1):
                red_factor = FLOAT_TYPE(1.0)
                time = FLOAT_TYPE(config.times[itime])
                tau = FLOAT_TYPE(config.end_jd - time)

                xc_ejec = FLOAT_TYPE(comet.xarr[itime])
                yc_ejec = FLOAT_TYPE(comet.yarr[itime])
                zc_ejec = FLOAT_TYPE(comet.zarr[itime])
                vxc_ejec = FLOAT_TYPE(comet.vxarr[itime])
                vyc_ejec = FLOAT_TYPE(comet.vyarr[itime])
                vzc_ejec = FLOAT_TYPE(comet.vzarr[itime])
                rc_ejec = FLOAT_TYPE(np.sqrt(xc_ejec ** 2 + yc_ejec ** 2 + zc_ejec ** 2))
                thetac_ejec = FLOAT_TYPE(comet.trueanarr[itime])

                if config.iejec_mode == 3:
                    subsolat = FLOAT_TYPE(np.arcsin(
                        np.sin(config.nuc_inc) * np.sin(config.nuc_phi + thetac_ejec)))
                    config.subsolar_file.write(
                        f"{time - config.per_jd:15.5e} {rc_ejec:15.5e} {thetac_ejec:15.5e} {subsolat / TORAD:15.5e}\n")

                xdtime = FLOAT_TYPE(time - config.per_jd)
                if comet.ec < 1.0:
                    while xdtime < -comet.half_per:
                        xdtime += comet.orb_per
                    while xdtime > comet.half_per:
                        xdtime -= comet.orb_per

                dmdtl, power, vfac, radmin, radmax = DustTail._interp5(xdtime, config)
                dmdt = FLOAT_TYPE(10.0 ** dmdtl)
                xmass = FLOAT_TYPE(dmdt * config.deltat * 86400.0)
                skip_particles = False

                if config.iejec_mode == 3 and config.isun == 1:
                    nshadow = 0
                    for kk in range(config.nevent):
                        random1 = FLOAT_TYPE(np.random.random())
                        area_longitude = FLOAT_TYPE(config.lon_min +
                                                    (config.lon_max - config.lon_min) * random1)
                        cl1 = FLOAT_TYPE(np.cos(HALFPI + config.lat_max))
                        cl2 = FLOAT_TYPE(np.cos(HALFPI + config.lat_min))
                        random2 = FLOAT_TYPE(np.random.random())
                        area_latitude = FLOAT_TYPE(np.arccos(
                            (cl2 - cl1) * random2 + cl1) - HALFPI)
                        _, ur, _, _ = DustTail._anisot_dir2(
                            thetac_ejec, time, area_latitude, area_longitude, config)
                        cosz = FLOAT_TYPE(-ur)
                        if cosz < 0.0:
                            nshadow += 1

                    if nshadow == config.nevent:
                        print(f"The entire prescribed area is in darkness at: "
                              f"{time - config.per_jd:10.3f} days to perihelion")
                        dmdt = FLOAT_TYPE(0.0)
                        xmass = FLOAT_TYPE(0.0)
                        skip_particles = True
                    else:
                        red_factor = FLOAT_TYPE(1.0 - float(nshadow) / float(config.nevent))
                        print(f"  Applying red_factor={red_factor:.4f} to dust production")
                        dmdt *= red_factor
                        xmass *= red_factor

                self.total_mass += xmass
                config.dustloss_file.write(f"    {time - config.per_jd:10.5E}     {rc_ejec:10.5E}     {dmdt:10.5E}\n")

                if skip_particles:
                    continue

                r1log = FLOAT_TYPE(np.log10(radmin))
                r2log = FLOAT_TYPE(np.log10(radmax))
                steprlog = FLOAT_TYPE((r2log - r1log) / float(config.nsizes))

                for irlog in range(config.nsizes):
                    rlog = FLOAT_TYPE(r1log + irlog * steprlog)
                    r1 = FLOAT_TYPE(10.0 ** rlog)
                    r2 = FLOAT_TYPE(10.0 ** (rlog + steprlog))
                    rad = FLOAT_TYPE(np.sqrt(r1 * r2))
                    eneint = FLOAT_TYPE(3.0 * xmass * powerint(r1, r2, power) /
                                        (FOURPI * config.pden * powerint(radmin, radmax, power + 3.0)))
                    beta = FLOAT_TYPE(1.191e-3 / (2.0 * config.pden * rad))
                    vej = FLOAT_TYPE(config.v0 * vfac * beta ** config.gamma * rc_ejec ** config.kappa)
                    xit, yit, zit = FLOAT_TYPE(0.0), FLOAT_TYPE(0.0), FLOAT_TYPE(0.0)

                    for k in range(config.nevent):
                        vxit, vyit, vzit = DustTail._get_ejection_velocity(
                            config.iejec_mode, vej, xc_ejec, yc_ejec, zc_ejec,
                            thetac_ejec, time, config
                        )
                        heliorbit = HelioOrbit()  # New instance for sequential execution
                        heliorbit.set_matrices(config.helio_matrix, None)
                        heliorbit.set_params_nm(
                            config.delta, config.nmpar, config.nmpar1, config.nmpar2,
                            config.nmpar3, config.nmpar4, config.nmpar5, config.nmpar6,
                            config.nmpar7, config.nmpar8
                        )
                        npar, mpar, lpar = heliorbit.heliorbit(
                            config.per_jd, config.tc, tau, QP, config.pden, rad,
                            xc_ejec, yc_ejec, zc_ejec, vxc_ejec, vyc_ejec, vzc_ejec,
                            config.rcobs, config.thetacobs, xit, yit, zit, vxit, vyit, vzit
                        )
                        psang_rad = FLOAT_TYPE(config.psang * TORAD)
                        cos_psang = FLOAT_TYPE(np.cos(psang_rad))
                        sin_psang = FLOAT_TYPE(np.sin(psang_rad))
                        nrot = FLOAT_TYPE(npar * cos_psang - mpar * sin_psang)
                        mrot = FLOAT_TYPE(npar * sin_psang + mpar * cos_psang)
                        npar, mpar = nrot, mrot

                        if config.igrapho == 1 and plot_handler is not None and plot_handler.available:
                            random1 = FLOAT_TYPE(np.random.random())
                            if random1 > config.pcp:
                                plot_handler.add_particle(npar, mpar)

                        if config.iprn == 1:
                            config.nm_file.write(
                                f"{npar:15.5e} {mpar:15.5e} {lpar:15.5e} {rad:15.5e} "
                                f"{itime:6d} {config.times[itime]:18.4f}\n")

                        if (npar >= config.nmin and npar <= config.nmax and
                                mpar >= config.mmin and mpar <= config.mmax):
                            ii = int((npar - config.nmin) * config.augx)
                            jj = int((mpar - config.mmin) * config.augy)
                            if 0 <= ii < config.nx and 0 <= jj < config.ny:
                                flux_out = FLOAT_TYPE(config.cte_part * rad ** 2)
                                particle_contribution = FLOAT_TYPE(flux_out * eneint / float(config.nevent))
                                config.flux[ii, jj] += particle_contribution
                                depth_contribution = FLOAT_TYPE(
                                    2.0 * eneint * PI * (rad / (config.grdsiz * 1.0e3)) ** 2 / float(config.nevent))
                                self.opt_depth[ii, jj] += depth_contribution
                                if lpar > 0.0 and ii == config.inuc and jj == config.jnuc:
                                    self.opt_depth_nuc += depth_contribution

    def apply_convolution(self, flux_array, sfwhm):
        """
        Apply Gaussian convolution to the image.

        Args:
            flux_array: Image array to convolve
            sfwhm: Full width at half maximum for Gaussian kernel
        """
        nx, ny = flux_array.shape
        sigma = FLOAT_TYPE(sfwhm / 2.35482)
        s2 = FLOAT_TYPE(sigma ** 2)
        ts = FLOAT_TYPE(10.0 * sigma)
        its = int(ts)
        itest = its % 2
        w = its + (1 if itest == 0 else 0)
        wbar = (w - 1) // 2
        garr = np.zeros(w, dtype=FLOAT_TYPE)
        sum_garr = FLOAT_TYPE(0.0)

        for i in range(w):
            x = FLOAT_TYPE(i - wbar)
            garr[i] = FLOAT_TYPE(np.exp(-x * x / (2.0 * s2)))
            sum_garr += garr[i]

        for i in range(w):
            garr[i] = FLOAT_TYPE(garr[i] / sum_garr)

        aux = np.zeros((nx, ny), dtype=FLOAT_TYPE)
        convolu = np.zeros((nx, ny), dtype=FLOAT_TYPE)

        for iy in range(ny):
            for ix in range(nx):
                val = FLOAT_TYPE(0.0)
                for i in range(w):
                    idx = ix - wbar + i
                    if 0 <= idx < nx:
                        val += FLOAT_TYPE(garr[i] * flux_array[idx, iy])
                aux[ix, iy] = val

        for iy in range(ny):
            for ix in range(nx):
                val = FLOAT_TYPE(0.0)
                for i in range(w):
                    idy = iy - wbar + i
                    if 0 <= idy < ny:
                        val += FLOAT_TYPE(garr[i] * aux[ix, idy])
                convolu[ix, iy] = val

        for ix in range(nx):
            for iy in range(ny):
                flux_array[ix, iy] = FLOAT_TYPE(convolu[ix, iy])

    def calculate_afrho_mag(self, flux, inuc, jnuc, grdsiz, cte_mag, magsun, scale, delta, rcobs, rho_ap=None):
        """
        Compute magnitude and Afrho of modeled tail within aperture.

        Args:
            flux: Flux array
            inuc: Nucleus x position
            jnuc: Nucleus y position
            grdsiz: Grid size in km/pixel
            cte_mag: Magnitude constant
            magsun: Sun magnitude
            scale: Image scale in arcsec/pixel
            delta: Distance to observer in AU
            rcobs: Heliocentric distance in AU
            rho_ap: Override aperture radius in km

        Returns:
            tuple: (afrho, afrho_0, mag) in meters and magnitude
        """
        rho_ap_km = rho_ap
        fcomet = FLOAT_TYPE(0.0)
        nx, ny = flux.shape

        for i in range(nx):
            for j in range(ny):
                radion = FLOAT_TYPE(
                    grdsiz * np.sqrt(float(inuc - i) ** 2 + float(jnuc - j) ** 2))
                if radion <= rho_ap_km:
                    fimag = flux[i, j]
                    if fimag <= 0.0:
                        fimag = FLOAT_TYPE(1.0e-25)
                    xmag = FLOAT_TYPE(cte_mag + magsun - 2.5 * np.log10(fimag) - 2.5 * np.log10(
                        scale * scale))
                    fcomet += FLOAT_TYPE(10 ** (-xmag / 2.5))

        if fcomet > 0.0:
            mag = FLOAT_TYPE(-2.5 * np.log10(fcomet))
        else:
            mag = FLOAT_TYPE(99.9)

        fsun = FLOAT_TYPE(10 ** (-magsun / 2.5))
        dist_factor = FLOAT_TYPE(2.0 * delta * AUKM * rcobs / rho_ap_km)
        cte = FLOAT_TYPE(dist_factor ** 2)
        flux_ratio = FLOAT_TYPE(fcomet / fsun)
        afrho = FLOAT_TYPE(cte * flux_ratio * rho_ap_km * 1.0e3)

        xp = FLOAT_TYPE(0.0)
        correc = FLOAT_TYPE(10 ** (SCHLEICHER_COEFFS[0] +
                                   xp * (SCHLEICHER_COEFFS[1] +
                                         xp * (SCHLEICHER_COEFFS[2] +
                                               xp * (SCHLEICHER_COEFFS[3] +
                                                     xp * (SCHLEICHER_COEFFS[4] +
                                                           xp * (SCHLEICHER_COEFFS[5] +
                                                                 xp * SCHLEICHER_COEFFS[6])))))))
        afrho_0 = FLOAT_TYPE(afrho / correc)

        return afrho, afrho_0, mag

    @staticmethod
    def _interp5(xdtime, config):
        """
        Simple linear interpolation on the dmdt_vel_power_rmin_rmax.dat file.

        Args:
            xdtime: Time value to interpolate
            config: Configuration object with dust loss rate data

        Returns:
            tuple: (dmdtl, power, vfac, radmin, radmax)
        """
        if xdtime < config.dtime[0] or xdtime > config.dtime[config.ninputs - 1]:
            print(f"Time to perihelion={xdtime} is out of bounds, using nearest value")
            if xdtime < config.dtime[0]:
                return (FLOAT_TYPE(config.dmdtlog[0]),
                        FLOAT_TYPE(config.powera[0]),
                        FLOAT_TYPE(config.velfac[0]),
                        FLOAT_TYPE(config.radiomin[0]),
                        FLOAT_TYPE(config.radiomax[0]))
            else:
                return (FLOAT_TYPE(config.dmdtlog[config.ninputs - 1]),
                        FLOAT_TYPE(config.powera[config.ninputs - 1]),
                        FLOAT_TYPE(config.velfac[config.ninputs - 1]),
                        FLOAT_TYPE(config.radiomin[config.ninputs - 1]),
                        FLOAT_TYPE(config.radiomax[config.ninputs - 1]))

        if xdtime == config.dtime[config.ninputs - 1]:
            return (FLOAT_TYPE(config.dmdtlog[config.ninputs - 1]),
                    FLOAT_TYPE(config.powera[config.ninputs - 1]),
                    FLOAT_TYPE(config.velfac[config.ninputs - 1]),
                    FLOAT_TYPE(config.radiomin[config.ninputs - 1]),
                    FLOAT_TYPE(config.radiomax[config.ninputs - 1]))

        left = 0
        right = config.ninputs - 1
        while right - left > 1:
            mid = (left + right) // 2
            if config.dtime[mid] > xdtime:
                right = mid
            else:
                left = mid

        i = right
        delta = FLOAT_TYPE(config.dtime[i] - config.dtime[i - 1])
        factor = FLOAT_TYPE((xdtime - config.dtime[i - 1]) / delta)

        dmdtl = FLOAT_TYPE(config.dmdtlog[i - 1] + factor * (config.dmdtlog[i] - config.dmdtlog[i - 1]))
        vfac = FLOAT_TYPE(config.velfac[i - 1] + factor * (config.velfac[i] - config.velfac[i - 1]))
        power = FLOAT_TYPE(config.powera[i - 1] + factor * (config.powera[i] - config.powera[i - 1]))
        radmin = FLOAT_TYPE(config.radiomin[i - 1] + factor * (config.radiomin[i] - config.radiomin[i - 1]))
        radmax = FLOAT_TYPE(config.radiomax[i - 1] + factor * (config.radiomax[i] - config.radiomax[i - 1]))

        return dmdtl, power, vfac, radmin, radmax

    @staticmethod
    def _anisot_dir2(thetac_ejec, time, area_latitude, area_longitude, config):
        """
        Calculate particle velocity components based on active area location.

        Args:
            thetac_ejec: True anomaly at ejection time
            time: Ejection time (JD)
            area_latitude: Latitude of active area (radians)
            area_longitude: Longitude of active area (radians)
            config: Configuration object

        Returns:
            tuple: (theta0, ur, utheta, uz) direction components
        """
        finu = FLOAT_TYPE(config.nuc_phi + thetac_ejec)
        cosfinu = FLOAT_TYPE(np.cos(finu))
        sinfinu = FLOAT_TYPE(np.sin(finu))
        cosi = FLOAT_TYPE(np.cos(config.nuc_inc))
        sini = FLOAT_TYPE(np.sin(config.nuc_inc))

        a11 = FLOAT_TYPE(cosfinu)
        a12 = FLOAT_TYPE(cosi * sinfinu)
        a13 = FLOAT_TYPE(sini * sinfinu)
        a21 = FLOAT_TYPE(-sinfinu)
        a22 = FLOAT_TYPE(cosi * cosfinu)
        a23 = FLOAT_TYPE(sini * cosfinu)
        a31 = FLOAT_TYPE(0.0)
        a32 = FLOAT_TYPE(sini)
        a33 = FLOAT_TYPE(-cosi)

        theta0 = FLOAT_TYPE(np.arctan(np.tan(thetac_ejec + config.nuc_phi) * cosi))
        theta0_t = FLOAT_TYPE(np.arctan(np.tan(config.nuc_phi) * cosi))
        theta = FLOAT_TYPE((TWOPI / config.period) * (time - config.per_jd) + theta0_t - theta0 + area_longitude)

        v1 = FLOAT_TYPE(np.cos(area_latitude) * np.cos(theta + theta0))
        v2 = FLOAT_TYPE(np.cos(area_latitude) * np.sin(theta + theta0))
        v3 = FLOAT_TYPE(np.sin(area_latitude))

        ur = FLOAT_TYPE(-(a11 * v1 + a12 * v2 + a13 * v3))
        utheta = FLOAT_TYPE(-(a21 * v1 + a22 * v2 + a23 * v3))
        uz = FLOAT_TYPE(-(a31 * v1 + a32 * v2 + a33 * v3))

        return theta0, ur, utheta, uz

    @staticmethod
    def _get_ejection_velocity(iejec_mode, vej, xc_ejec, yc_ejec, zc_ejec,
                               thetac_ejec, time, config):
        """
        Calculate ejection velocity based on the selected ejection mode.

        Args:
            iejec_mode: Ejection mode (1=isotropic, 2=sunward, 3=active areas)
            vej: Ejection velocity magnitude
            xc_ejec, yc_ejec, zc_ejec: Comet position
            thetac_ejec: True anomaly
            time: Ejection time
            config: Configuration object

        Returns:
            tuple: (vx, vy, vz) velocity components
        """
        if iejec_mode == 1:
            random1 = FLOAT_TYPE(np.random.random())
            random2 = FLOAT_TYPE(np.random.random())
            phi = FLOAT_TYPE(TWOPI * random1)
            theta = FLOAT_TYPE(np.arccos(2.0 * random2 - 1.0))
            vxit = FLOAT_TYPE(np.sin(theta) * np.cos(phi))
            vyit = FLOAT_TYPE(np.sin(theta) * np.sin(phi))
            vzit = FLOAT_TYPE(np.cos(theta))
            vxit = FLOAT_TYPE(CTEVEL * vej * vxit)
            vyit = FLOAT_TYPE(CTEVEL * vej * vyit)
            vzit = FLOAT_TYPE(CTEVEL * vej * vzit)

        elif iejec_mode == 2:
            while True:
                random1 = FLOAT_TYPE(np.random.random())
                random2 = FLOAT_TYPE(np.random.random())
                phi = FLOAT_TYPE(TWOPI * random1)
                theta = FLOAT_TYPE(np.arccos(2.0 * random2 - 1.0))
                vxit = FLOAT_TYPE(np.sin(theta) * np.cos(phi))
                vyit = FLOAT_TYPE(np.sin(theta) * np.sin(phi))
                vzit = FLOAT_TYPE(np.cos(theta))
                ejmod = FLOAT_TYPE(np.sqrt(xc_ejec ** 2 + yc_ejec ** 2 + zc_ejec ** 2))
                cosz = FLOAT_TYPE(-(vxit * xc_ejec + vyit * yc_ejec + vzit * zc_ejec) / ejmod)
                if cosz >= 0.0:
                    break
            vxit = FLOAT_TYPE(CTEVEL * vej * (cosz) ** config.expocos * vxit)
            vyit = FLOAT_TYPE(CTEVEL * vej * (cosz) ** config.expocos * vyit)
            vzit = FLOAT_TYPE(CTEVEL * vej * (cosz) ** config.expocos * vzit)

        else:
            while True:
                random1 = FLOAT_TYPE(np.random.random())
                area_longitude = FLOAT_TYPE(config.lon_min +
                                            (config.lon_max - config.lon_min) * random1)
                cl1 = FLOAT_TYPE(np.cos(HALFPI + config.lat_max))
                cl2 = FLOAT_TYPE(np.cos(HALFPI + config.lat_min))
                random2 = FLOAT_TYPE(np.random.random())
                area_latitude = FLOAT_TYPE(np.arccos((cl2 - cl1) * random2 + cl1) - HALFPI)
                _, ur, utheta, uz = DustTail._anisot_dir2(
                    thetac_ejec, time, area_latitude, area_longitude, config)
                cosz = FLOAT_TYPE(-ur)
                if not (cosz < 0.0 and config.isun == 1):
                    break
            vxop = FLOAT_TYPE(ur * np.cos(thetac_ejec) - utheta * np.sin(thetac_ejec))
            vyop = FLOAT_TYPE(ur * np.sin(thetac_ejec) + utheta * np.cos(thetac_ejec))
            vzop = FLOAT_TYPE(uz)
            from utils.coordinate_transforms import hpo_to_he
            vxit, vyit, vzit = hpo_to_he(vxop, vyop, vzop, config.helio_matrix)
            vxit = FLOAT_TYPE(CTEVEL * vej * (max(0.0, cosz)) ** config.expocos * vxit)
            vyit = FLOAT_TYPE(CTEVEL * vej * (max(0.0, cosz)) ** config.expocos * vyit)
            vzit = FLOAT_TYPE(CTEVEL * vej * (max(0.0, cosz)) ** config.expocos * vzit)

        return vxit, vyit, vzit

def dust_tail_worker(itime_start, itime_end, comet, config):
    """
    Worker function for parallel dust tail computation.

    Args:
        itime_start: Starting time step index
        itime_end: Ending time step index
        comet: Comet object
        config: SimulationConfig object

    Returns:
        dict: Local computation results
    """
    flux_local = np.zeros((config.nx, config.ny), dtype=FLOAT_TYPE)
    opt_depth_local = np.zeros((config.nx, config.ny), dtype=FLOAT_TYPE)
    total_mass_local = FLOAT_TYPE(0.0)
    opt_depth_nuc_local = FLOAT_TYPE(0.0)
    dustloss_data = []
    particles_to_plot = []
    subsolar_data = []

    heliorbit = HelioOrbit()
    heliorbit.set_matrices(config.helio_matrix, None)
    heliorbit.set_params_nm(
        config.delta, config.nmpar, config.nmpar1, config.nmpar2, config.nmpar3,
        config.nmpar4, config.nmpar5, config.nmpar6, config.nmpar7, config.nmpar8
    )

    if config.iprn == 1:
        nm_file = open(f'output/nm_part{os.getpid()}.dat', 'w')

    for itime in range(itime_start, itime_end - 1):
        red_factor = FLOAT_TYPE(1.0)
        time = FLOAT_TYPE(config.times[itime])
        tau = FLOAT_TYPE(config.end_jd - time)

        xc_ejec = FLOAT_TYPE(comet.xarr[itime])
        yc_ejec = FLOAT_TYPE(comet.yarr[itime])
        zc_ejec = FLOAT_TYPE(comet.zarr[itime])
        vxc_ejec = FLOAT_TYPE(comet.vxarr[itime])
        vyc_ejec = FLOAT_TYPE(comet.vyarr[itime])
        vzc_ejec = FLOAT_TYPE(comet.vzarr[itime])
        rc_ejec = FLOAT_TYPE(np.sqrt(xc_ejec ** 2 + yc_ejec ** 2 + zc_ejec ** 2))
        thetac_ejec = FLOAT_TYPE(comet.trueanarr[itime])

        if config.iejec_mode == 3:
            subsolat = FLOAT_TYPE(np.arcsin(
                np.sin(config.nuc_inc) * np.sin(config.nuc_phi + thetac_ejec)))
            subsolar_data.append((time - config.per_jd, rc_ejec, thetac_ejec, subsolat / TORAD))

        xdtime = FLOAT_TYPE(time - config.per_jd)
        if comet.ec < 1.0:
            while xdtime < -comet.half_per:
                xdtime += comet.orb_per
            while xdtime > comet.half_per:
                xdtime -= comet.orb_per

        dmdtl, power, vfac, radmin, radmax = DustTail._interp5(xdtime, config)
        dmdt = FLOAT_TYPE(10.0 ** dmdtl)
        xmass = FLOAT_TYPE(dmdt * config.deltat * 86400.0)
        skip_particles = False

        if config.iejec_mode == 3 and config.isun == 1:
            nshadow = 0
            for kk in range(config.nevent):
                random1 = FLOAT_TYPE(np.random.random())
                area_longitude = FLOAT_TYPE(config.lon_min +
                                            (config.lon_max - config.lon_min) * random1)
                cl1 = FLOAT_TYPE(np.cos(HALFPI + config.lat_max))
                cl2 = FLOAT_TYPE(np.cos(HALFPI + config.lat_min))
                random2 = FLOAT_TYPE(np.random.random())
                area_latitude = FLOAT_TYPE(np.arccos(
                    (cl2 - cl1) * random2 + cl1) - HALFPI)
                _, ur, _, _ = DustTail._anisot_dir2(
                    thetac_ejec, time, area_latitude, area_longitude, config)
                cosz = FLOAT_TYPE(-ur)
                if cosz < 0.0:
                    nshadow += 1

            if nshadow == config.nevent:
                print(f"The entire prescribed area is in darkness at: "
                      f"{time - config.per_jd:10.3f} days to perihelion")
                dmdt = FLOAT_TYPE(0.0)
                xmass = FLOAT_TYPE(0.0)
                skip_particles = True
            else:
                red_factor = FLOAT_TYPE(1.0 - float(nshadow) / float(config.nevent))
                print(f"  Applying red_factor={red_factor:.4f} to dust production")
                dmdt *= red_factor
                xmass *= red_factor

        total_mass_local += xmass
        dustloss_data.append((time - config.per_jd, rc_ejec, dmdt))

        if skip_particles:
            continue

        r1log = FLOAT_TYPE(np.log10(radmin))
        r2log = FLOAT_TYPE(np.log10(radmax))
        steprlog = FLOAT_TYPE((r2log - r1log) / float(config.nsizes))

        for irlog in range(config.nsizes):
            rlog = FLOAT_TYPE(r1log + irlog * steprlog)
            r1 = FLOAT_TYPE(10.0 ** rlog)
            r2 = FLOAT_TYPE(10.0 ** (rlog + steprlog))
            rad = FLOAT_TYPE(np.sqrt(r1 * r2))
            eneint = FLOAT_TYPE(3.0 * xmass * powerint(r1, r2, power) /
                                (FOURPI * config.pden * powerint(radmin, radmax, power + 3.0)))
            beta = FLOAT_TYPE(1.191e-3 / (2.0 * config.pden * rad))
            vej = FLOAT_TYPE(config.v0 * vfac * beta ** config.gamma * rc_ejec ** config.kappa)
            xit, yit, zit = FLOAT_TYPE(0.0), FLOAT_TYPE(0.0), FLOAT_TYPE(0.0)

            for k in range(config.nevent):
                vxit, vyit, vzit = DustTail._get_ejection_velocity(
                    config.iejec_mode, vej, xc_ejec, yc_ejec, zc_ejec,
                    thetac_ejec, time, config
                )
                npar, mpar, lpar = heliorbit.heliorbit(
                    config.per_jd, config.tc, tau, QP, config.pden, rad,
                    xc_ejec, yc_ejec, zc_ejec, vxc_ejec, vyc_ejec, vzc_ejec,
                    config.rcobs, config.thetacobs, xit, yit, zit, vxit, vyit, vzit
                )
                psang_rad = FLOAT_TYPE(config.psang * TORAD)
                cos_psang = FLOAT_TYPE(np.cos(psang_rad))
                sin_psang = FLOAT_TYPE(np.sin(psang_rad))
                nrot = FLOAT_TYPE(npar * cos_psang - mpar * sin_psang)
                mrot = FLOAT_TYPE(npar * sin_psang + mpar * cos_psang)
                npar, mpar = nrot, mrot

                if config.igrapho == 1:
                    random1 = FLOAT_TYPE(np.random.random())
                    if random1 > config.pcp:
                        particles_to_plot.append((npar, mpar))

                if config.iprn == 1:
                    nm_file.write(
                        f"{npar:15.5e} {mpar:15.5e} {lpar:15.5e} {rad:15.5e} "
                        f"{itime:6d} {config.times[itime]:18.4f}\n")

                if (npar >= config.nmin and npar <= config.nmax and
                        mpar >= config.mmin and mpar <= config.mmax):
                    ii = int((npar - config.nmin) * config.augx)
                    jj = int((mpar - config.mmin) * config.augy)
                    if 0 <= ii < config.nx and 0 <= jj < config.ny:
                        flux_out = FLOAT_TYPE(config.cte_part * rad ** 2)
                        particle_contribution = FLOAT_TYPE(flux_out * eneint / float(config.nevent))
                        flux_local[ii, jj] += particle_contribution
                        depth_contribution = FLOAT_TYPE(
                            2.0 * eneint * PI * (rad / (config.grdsiz * 1.0e3)) ** 2 / float(config.nevent))
                        opt_depth_local[ii, jj] += depth_contribution
                        if lpar > 0.0 and ii == config.inuc and jj == config.jnuc:
                            opt_depth_nuc_local += depth_contribution

    if config.iprn == 1:
        nm_file.close()

    return {
        'flux_local': flux_local,
        'opt_depth_local': opt_depth_local,
        'total_mass_local': total_mass_local,
        'opt_depth_nuc_local': opt_depth_nuc_local,
        'dustloss_data': dustloss_data,
        'particles_to_plot': particles_to_plot,
        'subsolar_data': subsolar_data if config.iejec_mode == 3 else []
    }