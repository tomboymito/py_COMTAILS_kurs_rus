"""
Version information for COMTAILS.

This module provides version and build information that is used throughout
the application to provide consistent version reporting.
"""


# Version information
VERSION = "1.0.0"
VERSION_DATE = "2025 May 05"
BUILD_INFO = f"COMTAILS - Comet Dust Tail Simulation v{VERSION} ({VERSION_DATE})"


def print_version_info():
    """
    Print version information to console.
    """
    print(f"\n{BUILD_INFO}")
    print("=" * len(BUILD_INFO))


def print_citation_info():
    """
    Print citation information to console.
    """
    citation = (" This code was ported from the FORTRAN serial version of   COMTAILS to Python by Rafael Morales and Nicolás Robles,\n of the Instituto de Astrofísica de  Andalucía. If you use this version of the code in your research,\n please make the appropriate acknowledgement to Rafael Morales and Nicolás Robles, and cite the source of the algorithm:\n\n"
        "Moreno, F. (2025). COMetary dust TAIL Simulator (COMTAILS): A computer code to generate comet dust tail \n brightness images. Astronomy and Astrophysics, Volume:695(2025), Article:A263."
    )

    print("\nCitation Information:")
    print("-" * 20)
    print(citation)
    print("-" * 20)
