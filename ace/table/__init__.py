# -*- coding: utf-8 -*-
"""Hash table module.

This module collects hash table generation methods.

Authors:
    Fangzhou Li: https://github.com/fangzhouli

"""

from .load import load_amino_acid_table
from .retrieve import get_fastas
from .generate import generate_table_entries

__all__ = [
    'get_fastas',
    'load_amino_acid_table',
    'generate_table_entries'
]
