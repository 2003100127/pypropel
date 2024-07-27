__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"


class Single:

    def __init__(
            self,
            msa,
    ):
        self.msa = msa
        self.msa_row = len(self.msa)

    def extract(
            self,
            i : int,
    ):
        """
        ith column in a MSA is extracted

        Examples
        --------
        output - ['A', 'A', 'A', ..., 'S']

        Parameters
        ----------
        i : int
            ith column in MSA

        Returns
        -------
        1d array
            row: 1; col: the number of multiple protein sequences.

        """
        return [row[i] for row in self.msa]

    def extract1(
            self,
            i : int,
    ):
        col_i = []
        for j in range(self.msa_row):
            col_i.append(self.msa[j][i])
        return col_i