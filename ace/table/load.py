# -*- coding: utf-8 -*-
"""Hash table generation function.

This module loads the hash table for users. For the first-time user, it
initializes an empty data/ht.pkl locally.

Attributes:
    PATH_DATA (str): The absolute path to the data folder.
    PATH_HT (str): The absolute path to the hash table pickle file.

Authors:
    Fangzhou Li: https://github.com/fangzhouli

"""

from os import path, mkdir
import pickle
import pandas as pd
from ace.table.retrieve import get_fastas
from ace.table.generate import generate_table_entries

PATH_DATA = path.abspath(path.dirname(__file__)) + '/data'
PATH_HT = path.abspath(path.dirname(__file__)) + '/data/ht.pkl'


def load_amino_acid_table(data):
    """Return a amino acid profile table that will be used to calculate AA
    composition. If some proteins are not in the hash table, the function will
    download them from UniProt. If certain proteins are not found on the
    website, it will generate a file to notify the user.

    Args:
        data (pd.DataFrame): An input proteomics data.

    Returns:
        aac_table (pd.DataFrame): Each row is a protein, and each column is an
            amino acid.

    """
    if not path.exists(PATH_DATA):
        mkdir(PATH_DATA)
        with open(PATH_HT, 'wb') as f:
            pickle.dump({}, f)

    proteins = data['Protein'].unique()
    accs_table = [protein.split('|')[1] for protein in proteins]
    accs_hashing = []
    ht = None
    aac_table = None

    with open(PATH_HT, 'rb') as f:
        ht = pickle.load(f)

    # Check the existence of new protein entries.
    for i, protein in enumerate(proteins):
        acc = protein.split('|')[1]
        if acc not in ht:
            accs_hashing.append((i, acc))

    # Generate new hash table entries if needed.
    if accs_hashing:
        fastas = get_fastas([x[1] for x in accs_hashing])
        accs_used = generate_table_entries(
            [x[1] for x in accs_hashing], fastas)

        # Substitute inconsistent accessions by used ones.
        for i, acc in enumerate(accs_hashing):
            idx, acc_old = acc
            accs_table[idx] = accs_used[i]

        with open(PATH_HT, 'rb') as f:
            ht = pickle.load(f)

    aac_table = pd.DataFrame(
        index=accs_table,
        columns=ht[next(iter(ht))].keys())
    aac_table = aac_table.apply(
        lambda row: pd.Series(ht[row.name]),
        axis=1)
    return aac_table
