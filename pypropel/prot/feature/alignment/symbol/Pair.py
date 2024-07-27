__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"


class Pair:

    def __init__(
            self,
            msa,
    ):
        self.msa = msa
        self.msa_row = len(self.msa)

    def extract(
            self,
            i : int,
            j : int,
    ):
        """
        The ith column and jth column in a MSA are extracted

        Examples
        --------
            output - (['A', 'A', 'A', ..., 'S'], ['T', 'T', 'Q', ..., 'T'])

        Parameters
        ----------
        i
            ith column in MSA
        j
            jth column in MSA

        Returns
        -------
        tuple
            two 1d array

        """
        col_i = []
        col_j = []
        for k in range(self.msa_row):
            col_i.append(self.msa[k][i])
            col_j.append(self.msa[k][j])
        return col_i, col_j

    def combine(
            self,
            i : int,
            j : int,
    ):
        """
        Pair consists of amino acid in ith column and amino acid
           in jth column in a MSA.

        Examples
        --------
            output - [('A', 'T'), ('-', 'A'), ..., ('-', 'T')]

        Parameters
        ----------
        i
            ith column in MSA
        j
            jth column in MSA

        Returns
        -------
        list : 1d array

        """
        col_i, col_j = self.extract(i, j)
        pairs = []
        for k in range(len(col_i)):
            pairs.append((col_i[k], col_j[k]))
        return pairs