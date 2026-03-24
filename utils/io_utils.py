"""Утилиты ввода-вывода и работы с файлами для COMTAILS."""
import os
import shutil

def reset_directory(dir_path):
    """Пересоздать каталог: удалить существующий и создать заново."""
    # Remove the directory if it exists
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
        print(f"Удалён существующий каталог: {dir_path}")

    # Create the directory again
    os.makedirs(dir_path)
    print(f"Создан каталог: {dir_path}")

