"""
Информация о версии COMTAILS.

Модуль предоставляет сведения о версии и сборке,
используемые для единообразного вывода в интерфейсе и CLI.
"""


# Version information
VERSION = "1.0.0"
VERSION_DATE = "2025 May 05"
BUILD_INFO = f"COMTAILS — моделирование пылевого хвоста кометы v{VERSION} ({VERSION_DATE})"


def print_version_info():
    """Вывести информацию о версии в консоль."""
    print(f"\n{BUILD_INFO}")
    print("=" * len(BUILD_INFO))


def print_citation_info():
    """Вывести рекомендации по научному цитированию."""
    citation = (
        " Этот код перенесён с последовательной версии COMTAILS на FORTRAN в Python\n"
        " Рафаэлем Моралесом и Николасом Роблесом (Instituto de Astrofísica de Andalucía).\n"
        " При использовании данной версии в научной работе, пожалуйста, укажите\n"
        " соответствующую благодарность авторам переноса и цитируйте источник алгоритма:\n\n"
        "Moreno, F. (2025). COMetary dust TAIL Simulator (COMTAILS): A computer code to generate comet dust tail\n"
        "brightness images. Astronomy and Astrophysics, Volume:695(2025), Article:A263."
    )

    print("\nИнформация для цитирования:")
    print("-" * 20)
    print(citation)
    print("-" * 20)
