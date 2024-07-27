__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np


class Binary:

    def __init__(self, msa):
        self.msa = msa
        self.msa_row = len(self.msa)
        self.msa_col = len(self.msa[0])

    def onehot(self):
        binary_matrix = [[0 for _ in range(self.msa_col * 21)] for _ in range(self.msa_row)]
        for i in range(self.msa_row):
            for j in range(self.msa_col):
                if 'A' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 0] = 1
                elif 'C' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 1] = 1
                elif 'D' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 2] = 1
                elif 'E' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 3] = 1
                elif 'F' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 4] = 1
                elif 'G' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 5] = 1
                elif 'H' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 6] = 1
                elif 'I' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 7] = 1
                elif 'K' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 8] = 1
                elif 'L' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 9] = 1
                elif 'M' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 10] = 1
                elif 'N' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 11] = 1
                elif 'P' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 12] = 1
                elif 'Q' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 13] = 1
                elif 'R' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 14] = 1
                elif 'S' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 15] = 1
                elif 'T' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 16] = 1
                elif 'V' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 17] = 1
                elif 'W' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 18] = 1
                elif 'Y' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 19] = 1
                elif '-' == self.msa[i][j]:
                    binary_matrix[i][j * 21 + 20] = 1
        binary_matrix = np.array(binary_matrix)
        return binary_matrix