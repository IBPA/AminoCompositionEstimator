# -*- coding: utf-8 -*-
"""FASTA retrieval functions.

This module retrieves FASTA files that are required to create the hash table.

Authors:
    Fangzhou Li - https://github.com/fangzhouli

"""

from http.client import HTTPSConnection
from multiprocessing import Pool
import logging


def _get_fasta_worker(acc):
    """A worker that retrieves the FASTA string of a given accession.

    Args:
        acc (str): The accession.

    Returns:
        res (str): The FASTE content.

    """
    conn = HTTPSConnection("www.uniprot.org")
    conn.request('GET', f"/uniprot/{acc}.fasta")
    res = conn.getresponse().read()
    if not res:
        conn.request('GET', f"/uniprot/?query={acc}&format=fasta")
        res = conn.getresponse().read()
    conn.close()
    return res.decode('utf-8')


def get_fastas(accs):
    """Retrieve FASTAs for given accessions from UniProt.

    Args:
        accs (list::str): A list of accessions.

    Returns:
        responses (list::str): A list of FASTA strings.

    """
    not_found = []

    with Pool(100) as p:
        responses = p.map(_get_fasta_worker, accs)

    for i, res in enumerate(responses):
        if not res:
            not_found.append(accs[i])

    if not_found:
        logging.info((
            "These protein entries' FASTAs were not able to be automatically "
            "retrieved. Please manually retrieve them."))
    return responses
