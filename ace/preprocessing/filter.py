# -*- coding: utf-8 -*-
"""Filtering functions.

This module includes filtering functions to exclude proteins with certain
patterns in the proteomics data.

Authors:
    Fangzhou Li: https://github.com/fangzhouli

"""

import logging


def filter_by_substr(data, substr):
    """Remove all proteins that include the specified substrings.

    Args:
        data (pd.DataFrame): Each row contains a peptide sequnce.

    Returns:
        (pd.DataFrame): Data with crap proteins removed.

    """
    logging.info("Removing filtered proteins from the data...")

    len_old = len(data)
    data = data[~data['Protein'].str.contains('|'.join(substr))]

    logging.debug("Removed {} proteins.".format(len_old - len(data)))
    return data
