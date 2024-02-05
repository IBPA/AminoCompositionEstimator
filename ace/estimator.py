# -*- coding: utf-8 -*-
"""

This module contains the API of estimating amino acid compositions of
IIFH food samples.

Authors:
    Fangzhou Li: https://github.com/fangzhouli

"""

import numpy as np
import pandas as pd
from ace import preprocessing
from ace.quant import topn
from ace.table import load_amino_acid_table


class Estimator:
    """The amino acid composition estimator using proteomics MS peak area data.

    Args:
        remove_zeroes (bool): If true, remove zero abundances in the proteomics
            data. Default is true.
        remove_duplicates (bool): If true, remove non-unique sequences occurred
            in different proteins. Default is true.
        merge_substr (bool): If true, merge all substring sequences in each
            protein.

    """

    def __init__(
            self,
            remove_zeroes=True,
            remove_duplicates=True,
            merge_substr=True):
        self.remove_zeroes = remove_zeroes
        self.remove_duplicates = remove_duplicates
        self.merge_substr = merge_substr

    def preprocess(self, data):
        """Preprocesses the proteomics data.

        Args:
            data: The peptide abundances.

        Returns:
            (pd.DataFrame): The preprocessed data.

        """
        if self.remove_zeroes:
            data = preprocessing.impute_zeroes_with_nan(data)
        if self.remove_duplicates:
            data = preprocessing.remove_non_unique_seq(data)
        if self.merge_substr:
            data = preprocessing.merge_seq_substr(data)
        return data

    def estimate(self, data):
        """Estimate the AA composition for proteomics data.

        Args:
            data (pd.DataFrame): The peptide abundances.

        Returns:
            aac (pd.DataFrame): The estimated amino acid composition.
            pai (pd.DataFrame): The estimated PAI.

        """
        data = self.preprocess(data)
        aac_table = load_amino_acid_table(data)
        pai = topn(data)

        # Normalize so that each column adds to 1.
        pai_norm = pai.div(pai.apply(np.sum, axis=0)).fillna(0)
        aac_table_norm = aac_table.div(aac_table.apply(np.sum, axis=1), axis=0)
        aac_table_norm.index = pai_norm.index
        aac = aac_table_norm.T.dot(pai_norm)
        return aac, pai


# if __name__ == '__main__':
#     filename = "~/Desktop/projects/ace/data/processed/tomato.csv"

#     est = Estimator()
#     data = pd.read_csv(filename)
#     aac, pai = est.estimate(data)
#     aac.to_csv('tomato.csv')
