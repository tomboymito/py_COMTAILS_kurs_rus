"""
I/O utilities for COMTAILS simulation.

This module provides input/output utility functions for file management.
"""
import os
import shutil

def reset_directory(dir_path):
    """
    Reset a directory by removing it and recreating it.

    Args:
        dir_path: Path to directory
    """
    # Remove the directory if it exists
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
        print(f"Deleted existing directory: {dir_path}")

    # Create the directory again
    os.makedirs(dir_path)
    print(f"Created directory: {dir_path}")


