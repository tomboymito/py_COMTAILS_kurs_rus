# COMTAILS — симулятор пылевого хвоста кометы

COMTAILS — Python-реализация модели формирования пылевого хвоста кометы. Проект основан на оригинальном FORTRAN-коде Fernando Moreno (IAA-CSIC) и адаптирован в объектно-ориентированной архитектуре с сохранением научной логики расчёта.

## Возможности

- Монте-Карло моделирование динамики пылевых частиц.
- Расчёт орбиты с использованием решения уравнения Кеплера.
- Эфемериды через JPL Horizons.
- Подмешивание звёздного поля (Gaia EDR3).
- Научные выходы в FITS (`tail_sdu.fits`, `tail_mag.fits`, `OPT_DEPTH.fits`).
- Расчёт фотометрических параметров `Afrho` и интегральной звёздной величины.
- Визуализация распределения частиц и fallback preview изображения.
- Два режима запуска: CLI и GUI.

## Установка

### Зависимости

- Python 3.8+
- `numpy`
- `astropy`
- `pygame`
- `requests`

Установка:

```bash
pip install numpy astropy pygame requests
```

## Запуск

### 1) CLI-режим (по умолчанию)

```bash
python main.py
```

С пользовательскими путями:

```bash
python main.py --input-dir input --output-dir output --config TAIL_INPUTS.dat --dust-profile dmdt_vel_power_rmin_rmax.dat
```

### 2) GUI-режим

```bash
python main.py --gui
```

GUI + автоматический запуск расчёта сразу после открытия окна:

```bash
python main.py --gui --gui-autorun
```

В GUI доступны:

- **Журнал вычислений** (пошаговый лог расчёта);
- **Итоговый блок результата** (`Afrho`, звёздная величина, масса пыли, путь к изображению);
- **Показ итогового изображения**.

## Логика итогового изображения и fallback

После расчёта программа ищет итоговое изображение в следующем порядке:

1. `output/dust_particles.png`
2. `output/tail_preview.ppm`
3. `output/tail_preview.png`

Если `dust_particles.png` не создан, формируется fallback-превью из расчётного потока яркости:

- обязательно `tail_preview.ppm`;
- дополнительно `tail_preview.png` (если доступно сохранение PNG через `pygame`).

## Выходные файлы

В каталоге `output` создаются:

- `tail_sdu.fits` — карта яркости пылевого хвоста (в единицах солнечного диска);
- `tail_mag.fits` — карта яркости в `mag/arcsec²`;
- `OPT_DEPTH.fits` — карта оптической толщины;
- `afrho.dat` — эволюция параметров `Afrho`;
- `dustlossrate.dat` — темп потери пыли;
- `dust_particles.png` — визуализация частиц (если включён графический вывод);
- `tail_preview.ppm`, `tail_preview.png` — fallback-превью, если основное PNG отсутствует.

## Структура проекта

- `main.py` — точка входа, поддержка CLI и `--gui`;
- `gui/app.py` — простой русскоязычный графический интерфейс;
- `simulation.py` — основной контроллер моделирования;
- `config.py` — чтение и хранение параметров;
- `models/` — физические модели кометы и пыли;
- `visualization/plot_handler.py` — графика частиц и подписи на русском;
- `visualization/preview_utils.py` — генерация fallback preview.

## Научная терминология

В пользовательской части используются астрофизически корректные термины:

- **Прямое восхождение (ПВ, RA)**;
- **Склонение (DEC)**;
- **перигелий / время до перигелия**;
- **оптическая толщина**, **звёздная величина**, **пылепроизводительность**.

## Лицензия

MIT, подробности см. в `LICENSE`.

## Цитирование

Если вы используете COMTAILS в научной работе, укажите происхождение модели и публикацию:

- Moreno, F. (2025). *COMetary dust TAIL Simulator (COMTAILS): A computer code to generate comet dust tail brightness images*. Astronomy and Astrophysics, 695, A263.
