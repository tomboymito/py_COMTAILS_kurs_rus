"""Модуль записи FITS-файлов для COMTAILS."""
import os
import numpy as np
import astropy.io.fits as fits

class FitsWriter:
    """Класс записи результатов моделирования в FITS с нужными заголовками."""

    def __init__(self):
        """Инициализировать модуль записи FITS."""
        # Ensure output directory exists
        os.makedirs('output', exist_ok=True)

    def write_fits_image(self, outimage, imagefit, object_name, start_jd, end_jd,
                         grdsiz, ntotmc, swap_axis_and_subtract_1_1=True):
        """Записать изображение в FITS.

        Параметры:
            outimage: путь к выходному FITS-файлу
            imagefit: двумерный массив изображения
            object_name: имя объекта (кометы)
            start_jd: начальная юлианская дата
            end_jd: дата наблюдения (юлианская)
            grdsiz: размер пикселя сетки, км/пикс
            ntotmc: общее число Монте-Карло событий
            swap_axis_and_subtract_1_1: поменять оси местами и сдвинуть матрицу на (-1, -1)
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
        print(f"Записан FITS-файл: {outimage}")
