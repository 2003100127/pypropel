__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np
from collections import Counter
from pypropel.prot.feature.alignment.symbol.Pair import Pair as pairalignsymbol
from pypropel.prot.sequence.Symbol import Symbol
from pypropel.util.ComputLib import ComputLib


class Pair:

    def __init__(self, msa):
        self.msa = msa
        self.computlib = ComputLib()
        self.msa_row = len(self.msa)
        self.pasp = pairalignsymbol(self.msa)
        self.aa_alphabet = Symbol().single(gap=True)

    def frequency_deprecated(self, i, j):
        """

        Parameters
        ----------
        i
            ith column in MSA
        j
            jth column in MSA

        Returns
        -------

        """
        aa_21x21 = self.computlib.combo2x2(self.aa_alphabet)
        res_pairs = self.pasp.combine(i, j)
        freq = {}
        for _, aa in enumerate(aa_21x21):
            freq[aa[0] + aa[1]] = 0
        for _, aa in enumerate(aa_21x21):
            for j, rr in enumerate(res_pairs):
                if np.array_equal(aa, rr):
                    freq[aa[0] + aa[1]] = freq[aa[0] + aa[1]] + 1
        for _, aa in enumerate(aa_21x21):
            freq[aa[0] + aa[1]] = freq[aa[0] + aa[1]] / self.msa_row
        return freq

    def frequency(self, i, j):
        """
        Fast frequency calculation.
        It provides methods to calculate pair frequency of two amino acids based on MSA.

        Parameters
        ----------
        i
            ith column in MSA
        j
            jth column in MSA

        Returns
        -------
            1d dict : dict

        """
        pairs = self.pasp.combine(i, j)
        aa_21x21 = self.computlib.combo2x2(self.aa_alphabet)
        count = {}
        for _, aa in enumerate(aa_21x21):
            count[aa[0] + aa[1]] = 0
        count_ = Counter(pairs).most_common()
        # print(count_)
        for _, aa in enumerate(count_):
            count[aa[0][0] + aa[0][1]] = aa[1]
        # print(count)
        for _, aa in enumerate(aa_21x21):
            count[aa[0] + aa[1]] = count[aa[0] + aa[1]] / self.msa_row
        freq = count
        return freq