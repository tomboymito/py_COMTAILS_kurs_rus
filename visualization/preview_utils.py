"""Утилиты для создания fallback-превью итогового изображения."""
from __future__ import annotations

import os
import numpy as np

try:
    import pygame as pg
except Exception:  # pygame опционален
    pg = None


def _normalize_to_uint8(image: np.ndarray) -> np.ndarray:
    """Нормализовать массив яркости в диапазон 0..255."""
    data = np.asarray(image, dtype=float)
    finite_mask = np.isfinite(data)

    if not np.any(finite_mask):
        return np.zeros_like(data, dtype=np.uint8)

    min_val = np.nanmin(data)
    max_val = np.nanmax(data)

    if max_val <= min_val:
        return np.zeros_like(data, dtype=np.uint8)

    scaled = (data - min_val) / (max_val - min_val)
    scaled = np.clip(scaled, 0.0, 1.0)
    return (scaled * 255.0).astype(np.uint8)


def save_flux_ppm(image: np.ndarray, output_path: str) -> str:
    """Сохранить изображение яркости в формате PPM (P6)."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    gray = _normalize_to_uint8(image)
    rgb = np.repeat(gray[:, :, np.newaxis], 3, axis=2)

    header = f"P6\n{rgb.shape[1]} {rgb.shape[0]}\n255\n".encode("ascii")
    with open(output_path, "wb") as f:
        f.write(header)
        f.write(rgb.tobytes())

    return output_path


def save_flux_png(image: np.ndarray, output_path: str) -> str | None:
    """Сохранить изображение яркости в PNG через pygame (если доступно)."""
    if pg is None:
        return None

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    gray = _normalize_to_uint8(image)
    rgb = np.repeat(gray[:, :, np.newaxis], 3, axis=2)

    try:
        if not pg.get_init():
            pg.init()
        surface = pg.surfarray.make_surface(np.transpose(rgb, (1, 0, 2)))
        pg.image.save(surface, output_path)
        return output_path
    except Exception:
        return None
