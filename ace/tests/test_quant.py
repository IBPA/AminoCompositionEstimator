import unittest
import numpy as np
import pandas as pd
from ace.quant import topn

class TestQuantMethods(unittest.TestCase):

    def test_top3(self):
        data = pd.DataFrame(
            {'Protein': ['a', 'a', 'a', 'a'],
             'Sequence': ['I', 'dont', 'really', 'care'],
             'col1': [1, 2, 3, 4],
             'col2': [4, 5, np.nan, np.nan],
             'col3': [np.nan, np.nan, np.nan, np.nan]})
        res = topn(data, 3)
        self.assertEqual(res['col1']['a'], 3)
        self.assertEqual(res['col2']['a'], 4.5)
        self.assertTrue(np.isnan(res['col3']['a']))


if __name__ == '__main__':
    unittest.main()
