# -*- coding: utf-8 -*-
"""Imputation functions.

This module includes imputation approaches for the proteomics data.

Authors:
    Fangzhou Li: https://github.com/fangzhouli

"""

import logging
import numpy as np


def impute_zeroes_with_nan(data):
    """Remove all zero-valued abundances in the proteomics data. This avoids
    Top3 mechanism to consider zeroes while calculating the mean.

    Args:
        data (pd.DataFrame): Each row contains a peptide sequnce.

    Returns:
        (pd.DataFrame): Data with zero-valued abundances removed.

    """
    logging.info("Imputing zero-valued abundances with NaN...")

    return data.replace({0: np.nan})
