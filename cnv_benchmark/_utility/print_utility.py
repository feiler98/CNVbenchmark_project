"""
------------------------------------------------------------------------------------------------------------------------

PRINT_UTILITY
-------------

Utility for the multiomics and facs data.
General functions which optimize handling of the data

------------------------------------------------------------------------------------------------------------------------
"""

# imports
# ______________________________________________________________________________________________________________________
import ast
import re
from pathlib import Path
from Bio import Align
import pyomics
# ______________________________________________________________________________________________________________________

def print_list(list_input: list, header_str: str):
    """
    Parameters
    ----------
    list_input: list
    header_str: str
        The text added above the displayed data
    """
    # box in the information for terminal display
    print("="*(14+len(header_str)))
    print(f"    // {header_str}")
    print("-" * (14 + len(header_str)))

    list_count_cells = []
    dict_group = {}
    for items in list_input:
        print(f"> {items}")

    print("=" * (14 + len(header_str)))