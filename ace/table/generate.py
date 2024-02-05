# -*- coding: utf-8 -*-
"""Hash table generation function.

This module creates a hash table for users. For the first-time user, it
generates a data/ht.pkl locally.

Attributes:
    AMINO_ACID_SYMBOLS (str): Amino acid sequence characters.

Authors:
    Fangzhou Li: https://github.com/fangzhouli

"""

from os import path
import pickle
import logging

AMINO_ACID_SYMBOLS = 'ABCDEFGHIKLMNPQRSTUVWXYZ'


def parse_fasta(fasta):
    """Parse a FASTA file and get accession and sequence.

    Args:
        fasta (str): A FASTA file in string format.

    Returns:
        protein_accession (str): A unique ID for the FASTA.
        sequence (str): A one-line sequence for the protein.

    """
    lines = fasta.split('\n')
    description = lines[0]
    if description[:3] not in ['>sp', '>tr']:
        logging.critical((
            f"Protein {description} is not one of the supported "
            f"formats."))
        raise ValueError("Only support '>sp' and '>tr'.")
    protein_accession = description.split('|')[1]
    sequence = ''.join(lines[1:])
    return (protein_accession, sequence)


def get_amino_acid_counts(protein_accession, seq):
    """Extract AA counts for a given protein sequence.

    Args:
        seq (str): A sequence consisted by only AA symbols.

    Returns:
        counts (dict): A dictionary with 20 AA symbol keys, and value of
            non-negative integer.

    """
    counts = {}
    for aa in AMINO_ACID_SYMBOLS:
        counts[aa] = 0

    for aa in seq:
        if aa not in AMINO_ACID_SYMBOLS:
            logging.critical(
                f"Protein {protein_accession} contains non-AA characters.")
        counts[aa] += 1
    return counts


def generate_table_entries(accs, fastas):
    """Generate a hash table with key of protein accession and value of amino
    acid counts.

    Args:
        fastas (list of str): A list of FASTA files.

    Returns:
        accs_used (list of str): A list of accessions that are consistent with
            UniProt.

    Generates:
        inconsistence.txt (file): A file that includes a list of accessions
            that search engine return but with different accessions provided by
            the user.

    """
    logging.info(f"{len(fastas)} new proteins found, hashing new entries.")

    ht = None
    accs_used = []
    accs_check = []
    with open(
            path.abspath(path.dirname(__file__)) + '/data/ht.pkl',
            'rb') as f:
        ht = pickle.load(f)

    for acc, fasta in zip(accs, fastas):
        protein_accession, sequence = parse_fasta(fasta)
        accs_used.append(protein_accession)
        if acc != protein_accession:
            accs_check.append((acc, protein_accession))
        ht[protein_accession] = get_amino_acid_counts(
            protein_accession, sequence)

    with open(
            path.abspath(path.dirname(__file__)) + '/data/ht.pkl', 'wb') as f:
        pickle.dump(ht, f)

    if accs_check:
        logging.warning(
            (f"{len(accs_check)} inconsistent accessions detected."
             f" See them in `inconsistence.txt`."))
        with open('inconsistence.txt', 'w') as f:
            f.write(str(accs_check))
    return accs_used
