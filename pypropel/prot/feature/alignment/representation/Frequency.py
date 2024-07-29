__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np
from collections import Counter
from pypropel.prot.feature.alignment.symbol.Single import Single as sglalignsymbol


class Frequency:

    def __init__(self, msa):
        self.msa = msa
        self.msa_col = len(self.msa[0])
        self.passingle = sglalignsymbol(self.msa)

    def calc_sgl_col(self, x):
        bases = self.passingle.extract(x)
        # print(bases[1])
        num_total = len(bases)
        freq_dictionary = Counter(bases)
        freq_all_bases = freq_dictionary.most_common()
        # print(freq_all_bases)
        template_num = len(freq_all_bases)
        # print(template_num)
        freq_single = [0 for _ in range(num_total)]
        # print(num_total)
        for i in range(template_num):
            for j in range(num_total):
                if np.array_equal(freq_all_bases[i][0], bases[j]):
                    freq_single[j] = round(freq_all_bases[i][1] / num_total, 4)
                    if bases[j] == '-':
                        freq_single[j] = 0
        return freq_single

    def matrix(self ):
        base_freq_matrix_T = []
        for i in range(self.msa_col):
            base_freq_matrix_T.append(self.calc_sgl_col(i))
        base_freq_matrix = np.transpose(np.array(base_freq_matrix_T))
        return base_freq_matrix


if __name__ == "__main__":
    p = Frequency(
        msa=msa
    )
    # print(p.calc_sgl_col(1))
    print(p.matrix())