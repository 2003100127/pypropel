__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import time
import numpy as np
import pandas as pd
from pypropel.prot.feature.alignment.frequency.Single import Single
from pypropel.prot.feature.alignment.frequency.Pair import Pair


class InformationTheory:
    
    def __init__(
            self,
            msa,
    ):
        self.msa = msa
        self.sfreq = Single(self.msa)
        self.pfreq = Pair(self.msa)
    
    def entropy(self, ):
        ent = dict()
        freq_column, _ = self.sfreq.columns()
        len_seq = len(self.msa[0])
        for i in range(len_seq):
            tmp = 0
            for aa in self.sfreq.aa_alphabet:
                if freq_column[aa][i] != 0:
                    tmp = tmp + freq_column[aa][i] * np.log10(freq_column[aa][i])
            ent[i + 1] = -tmp * 2
        return ent
    
    def entropy_gap(
            self, gap_thres=1
    ):
        ent = dict()
        freq_column, gap_compos = self.sfreq.columns()
        gap_compos_df = pd.DataFrame(gap_compos)
        filtered_indices = gap_compos_df.loc[gap_compos_df[0] < gap_thres].index
        # print(filtered_indices)
        for i in filtered_indices:
            tmp = 0
            for aa in self.sfreq.aa_alphabet:
                if freq_column[aa][i] != 0:
                    tmp -= freq_column[aa][i] * np.log2(freq_column[aa][i])
            ent[i + 1] = tmp
        return ent

    def entropy_(
            self,
            ent_dict,
            list_2d,
            window_aa_ids,
    ):
        start_time = time.time()
        list_2d_ = list_2d
        # ent_dict = self.entropy()
        for i, aa_win_ids in enumerate(window_aa_ids):
            for j in aa_win_ids:
                if j is None:
                    list_2d_[i].append(0)
                else:
                    list_2d_[i].append(ent_dict[j])
        # print(len(list_2d_[0]))
        # for i in list_2d_:
        #     if len(i) != 1:
        #         print(len(i))
        end_time = time.time()
        print('------> it entropy {time}s.'.format(time=end_time - start_time))
        return list_2d_

    def mi(self, i, j):
        """
        if freq_pair[aa1 + aa2] != 0, f1 != 0 and f2 != 0.

        Parameters
        ----------
        i
            ith column in MSA
        j
            jth column in MSA

        Returns
        -------

        """
        freq_column, _ = self.sfreq.columns()
        freq_pair = self.pfreq.frequency(i, j)
        res = 0
        for _, aa1 in enumerate(self.sfreq.aa_alphabet):
            f1 = freq_column[aa1][i]
            for _, aa2 in enumerate(self.sfreq.aa_alphabet):
                f2 = freq_column[aa2][j]
                if freq_pair[aa1 + aa2] != 0:
                    part = freq_pair[aa1 + aa2] / (f1 * f2)
                    # print(part)
                    res = res + freq_pair[aa1 + aa2] * np.log(part)
        return res