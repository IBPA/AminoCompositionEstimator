# -*- coding: utf-8 -*-
"""Preprocessing peptide sequences.

This module provides different proprocessing functions for peptides of the
proteomics data.

Authors:
    Fangzhou Li: https://github.com/fangzhouli

"""

from collections import defaultdict
import logging
import numpy as np
import pandas as pd


def remove_non_unique_seq(data):
    """Remove all peptides that are contained by different proteins from
    proteomics data.

    Args:
        data (pd.DataFrame): Each row contains a peptide sequnce.

    Returns:
        (pd.DataFrame): Data with each peptide only unique to one protein.

    """
    logging.info("Removing every non-unique sequence from the data...")
    # For each unique peptide sequence, count its occurrences in all proteins.
    seq_count = defaultdict(int)
    for seq in data['Sequence'].tolist():
        seq_count[seq] += 1
    # Ignore peptides that have occurances more than 1.
    total_seq_non_unique = 0
    seq_non_unique = []
    for seq, count in seq_count.items():
        if count > 1:
            total_seq_non_unique += count
            seq_non_unique.append(seq)
    logging.debug(
        (f"There are {len(seq_non_unique)} kinds of non-unique sequences. "
         f"Removed {total_seq_non_unique} sequences from the data."))
    return data[~data['Sequence'].isin(seq_non_unique)]


def merge_seq_substr(data):
    """Merge peptides which are substrings of another peptide.

    Args:
        data (pd.DataFrame): Peptides of a protein.

    Returns:
        (pd.DataFrame): Merged data.

    """
    logging.info("Merging non-unique peptides...")

    data_merged_list = []
    for protein in data['Protein'].unique().tolist():
        data_protein = data[data['Protein'] == protein]
        _data_protein = data_protein.copy()
        _data_protein.index = _data_protein['Sequence']
        _data_protein = _data_protein.iloc[:, 2:]
        # Each group is the longest peptide sequence of that group.
        #   For each peptide, check if it is a substring of any group. If
        #   yes, add this peptide to that group. If not, create a new group.
        seq_group = defaultdict(list)
        for seq in sorted(_data_protein.index.tolist(), key=len,
                          reverse=True):
            is_substr = False
            for key in seq_group.keys():
                if seq in key:
                    seq_group[key].append(seq)
                    is_substr = True
            if not is_substr:
                seq_group[seq].append(seq)
        # Generating a merged data.
        data_merged = data_protein[data_protein['Sequence'].isin(
            seq_group.keys())].copy()
        data_merged.iloc[:, 2:] = data_merged.apply(
            lambda row: _data_protein.loc[seq_group[row['Sequence']]].sum(
                axis=0),
            axis=1)
        data_merged_list.append(data_merged.replace(0, np.nan))
    return pd.concat(data_merged_list, axis='index')
