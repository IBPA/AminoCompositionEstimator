# -*- coding: utf-8 -*-
"""Topn quantification method.

This module implements Topn protein quantitation algorithm.

Authors:
    Fangzhou Li: https://github.com/fangzhouli

Reference:
    - Silva, Jeffrey C., et al. "Absolute quantification of proteins by LCMSE:
        a virtue of parallel MS acquisition." Molecular & Cellular Proteomics
        5.1 (2006): 144-156.

"""

import numpy as np
import pandas as pd


def topn_protein(data_protein, n=3):
    """Estimate PAI for a protein in each sample.

    Args:
        data_protein (pd.DataFrame): The processed protein data.
        n (int): Estimate PAI with N most signified peptides. By default n
            is 3.

    Returns:
        (pd.DataFrame): The PAI for the protein.

    """
    samples = list(data_protein.columns[2:])
    pai_protein = pd.Series(index=samples, dtype=np.float64)

    for sample in samples:
        data_sample = data_protein[sample].sort_values(ascending=False)
        data_sample_clean = data_sample[:n].dropna()
        if len(data_sample_clean):
            pai_protein[sample] = np.sum(data_sample_clean) / \
                len(data_sample_clean)
        else:
            pai_protein[sample] = np.nan
    return pai_protein


def topn(data, n=3):
    """Estimate the protein abundance index (PAI) for each protein in the each
    sample of the proteomics data.

    Args:
        data (pd.DataFrame): The MS peak areas proteomics data.
        n (int): Estimate PAI with N most signified peptides. By default n
            is 3.

    Returns:
        (pd.DataFrame): The quantified PAI. Each row is a protein, and each
            column is a sample.

    """
    proteins = data['Protein'].unique()
    pai = pd.DataFrame(index=proteins, columns=list(data.columns[2:]))

    for protein in proteins:
        pai.loc[protein] = topn_protein(data[data['Protein'] == protein], n)
    return pai
