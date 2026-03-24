"""Точка входа в COMTAILS.

Модуль запускает моделирование пылевого хвоста кометы
в консольном режиме или в графическом интерфейсе.
"""
import os
import sys
import argparse

from simulation import SimulationController
from utils.version import print_version_info, print_citation_info
from utils.io_utils import reset_directory


def parse_arguments():
    """Разобрать аргументы командной строки."""
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='COMTAILS: симулятор пылевого хвоста кометы')

    parser.add_argument('--input-dir', type=str, default='input',
                        help='Каталог со входными файлами')

    parser.add_argument('--output-dir', type=str, default='output',
                        help='Каталог для выходных файлов')

    parser.add_argument('--config', type=str, default='TAIL_INPUTS.dat',
                        help='Имя основного конфигурационного файла')

    parser.add_argument('--dust-profile', type=str, default='dmdt_vel_power_rmin_rmax.dat',
                        help='Имя файла профиля пылепотерь')

    parser.add_argument('--gui', action='store_true',
                        help='Запустить простой графический интерфейс (GUI)')
    parser.add_argument('--gui-autorun', action='store_true',
                        help='В режиме GUI сразу запустить расчёт после открытия окна')

    parser.add_argument('--validate', action='store_true',
                        help='Проверить результаты относительно эталонных значений')

    parser.add_argument('--expected-afrho', type=float, default=10.5,
                        help='Эталонное значение Afrho для проверки')

    parser.add_argument('--expected-mag', type=float, default=8.07,
                        help='Эталонная звёздная величина для проверки')

    parser.add_argument('--tolerance', type=float, default=0.1,
                        help='Допуск проверки (относительное отклонение)')

    return parser.parse_args()


def main():
    """Запустить моделирование COMTAILS."""
    # Parse command line arguments
    args = parse_arguments()

    # Run GUI mode if requested
    if args.gui:
        from gui import run_gui
        run_gui(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            config_file=args.config,
            dust_profile=args.dust_profile,
            auto_run=args.gui_autorun
        )
        run_gui()
        return

    # Check for required input files
    input_dir = args.input_dir
    required_files = [
        os.path.join(input_dir, args.config),
        os.path.join(input_dir, args.dust_profile)
    ]

    for file in required_files:
        if not os.path.exists(file):
            print(f"Ошибка: обязательный входной файл '{file}' не найден.")
            sys.exit(1)

    # Print version information
    print_version_info()

    # Create output directory
    reset_directory(args.output_dir)

    # Create and run simulation
    simulation = SimulationController()
    simulation.run(required_files)

    # Validate results if requested
    if args.validate:
        validation_result = SimulationController.validate_results(
            os.path.join(args.output_dir, "afrho.dat"),
            args.expected_afrho,
            args.expected_mag,
            args.tolerance
        )

        if not validation_result:
            print("Проверка не пройдена!")
            sys.exit(1)

    # Print version information again
    print_version_info()
    print_citation_info()


if __name__ == "__main__":
    main()
