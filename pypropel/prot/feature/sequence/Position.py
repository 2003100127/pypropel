__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from typing import List

import time
import numpy as np


class Position:

    def __init__(self, ):
        pass

    def deepconpred(self, ):
        """
        Sequence separation is encoded into binary representation with length of 9.

        References
        ----------
        Xiong, D., Zeng, J., & Gong, H. (2017). A deep learning framework
        for improving long-range residue–residue contact prediction using a
        hierarchical strategy. Bioinformatics, 33(17), 2675-2683.

        Returns
        -------
        1d dictionary : dict

        """
        num_sort = 10
        one_hot_binary = np.zeros((num_sort, num_sort))
        one_hot_binary[np.arange(num_sort), np.flipud(np.arange(num_sort))] = 1
        # print(one_hot_binary)
        encode_dict = {
            '24<=v<=28': one_hot_binary[0].tolist(),
            '29<=v<=33': one_hot_binary[1].tolist(),
            '34<=v<=38': one_hot_binary[2].tolist(),
            '39<=v<=43': one_hot_binary[3].tolist(),
            '44<=v<=48': one_hot_binary[4].tolist(),
            '49<=v<=58': one_hot_binary[5].tolist(),
            '59<=v<=68': one_hot_binary[6].tolist(),
            '69<=v<=78': one_hot_binary[7].tolist(),
            '79<=v<=88': one_hot_binary[8].tolist(),
            'v>=89': one_hot_binary[9].tolist(),
            'none': [0] * num_sort
        }
        return encode_dict

    def deepconpred_(self, list_2d, window_aa_ids):
        start_time = time.time()
        list_2d_ = list_2d
        # window_aa_ids_ = [i[0] + i[1] for i in window_aa_ids]
        encode_dict = self.deepconpred()
        for i, aa_win_ids in enumerate(window_aa_ids):
            res1 = aa_win_ids[0][0]
            res2 = aa_win_ids[1][0]
            if res1 == None or res2 == None:
                list_2d_[i] = list_2d_[i] + encode_dict['none']
            elif 24 <= abs(res1 - res2) <= 28:
                list_2d_[i] = list_2d_[i] + encode_dict['24<=v<=28']
            elif 29 <= abs(res1 - res2) <= 33:
                list_2d_[i] = list_2d_[i] + encode_dict['29<=v<=33']
            elif 34 <= abs(res1 - res2) <= 38:
                list_2d_[i] = list_2d_[i] + encode_dict['34<=v<=38']
            elif 39 <= abs(res1 - res2) <= 43:
                list_2d_[i] = list_2d_[i] + encode_dict['39<=v<=43']
            elif 44 <= abs(res1 - res2) <= 48:
                list_2d_[i] = list_2d_[i] + encode_dict['44<=v<=48']
            elif 49 <= abs(res1 - res2) <= 58:
                list_2d_[i] = list_2d_[i] + encode_dict['49<=v<=58']
            elif 59 <= abs(res1 - res2) <= 68:
                list_2d_[i] = list_2d_[i] + encode_dict['59<=v<=68']
            elif 69 <= abs(res1 - res2) <= 78:
                list_2d_[i] = list_2d_[i] + encode_dict['69<=v<=78']
            elif 89 <= abs(res1 - res2) <= 88:
                list_2d_[i] = list_2d_[i] + encode_dict['79<=v<=88']
            else:
                list_2d_[i] = list_2d_[i] + encode_dict['v>=89']
        end_time = time.time()
        print('------> absolute position: {time}s.'.format(time=end_time - start_time))
        return list_2d_

    def metapsicov(self):
        num_sort = 16
        one_hot_binary = np.zeros((num_sort, num_sort))
        one_hot_binary[np.arange(num_sort), np.flipud(np.arange(num_sort))] = 1
        # print(one_hot_cls)
        encode_dict = {
            'v<5': one_hot_binary[0].tolist(),
            'v=5': one_hot_binary[1].tolist(),
            'v=6': one_hot_binary[2].tolist(),
            'v=7': one_hot_binary[3].tolist(),
            'v=8': one_hot_binary[4].tolist(),
            'v=9': one_hot_binary[5].tolist(),
            'v=10': one_hot_binary[6].tolist(),
            'v=11': one_hot_binary[7].tolist(),
            'v=12': one_hot_binary[8].tolist(),
            'v=13': one_hot_binary[9].tolist(),
            '14<=v<18': one_hot_binary[10].tolist(),
            '18<=v<23': one_hot_binary[11].tolist(),
            '23<=v<28': one_hot_binary[12].tolist(),
            '28<=v<38': one_hot_binary[13].tolist(),
            '38<=v<48': one_hot_binary[14].tolist(),
            'v>=48': one_hot_binary[15].tolist(),
            'none': [0] * num_sort
        }
        return encode_dict

    def metapsicov_(self, list_2d, window_aa_ids):
        start_time = time.time()
        list_2d_ = list_2d
        # window_aa_ids_ = [i[0] + i[1] for i in window_aa_ids]
        encode_dict = self.metapsicov()
        for i, aa_win_ids in enumerate(window_aa_ids):
            res1 = aa_win_ids[0][0]
            res2 = aa_win_ids[1][0]
            if res1 == None or res2 == None:
                list_2d_[i] = list_2d_[i] + encode_dict['none']
            elif abs(res1 - res2) < 5:
                list_2d_[i] = list_2d_[i] + encode_dict['v<5']
            elif abs(res1 - res2) == 5:
                list_2d_[i] = list_2d_[i] + encode_dict['v=5']
            elif abs(res1 - res2) == 6:
                list_2d_[i] = list_2d_[i] + encode_dict['v=6']
            elif abs(res1 - res2) == 7:
                list_2d_[i] = list_2d_[i] + encode_dict['v=7']
            elif abs(res1 - res2) == 8:
                list_2d_[i] = list_2d_[i] + encode_dict['v=8']
            elif abs(res1 - res2) == 9:
                list_2d_[i] = list_2d_[i] + encode_dict['v=9']
            elif abs(res1 - res2) == 10:
                list_2d_[i] = list_2d_[i] + encode_dict['v=10']
            elif abs(res1 - res2) == 11:
                list_2d_[i] = list_2d_[i] + encode_dict['v=11']
            elif abs(res1 - res2) == 12:
                list_2d_[i] = list_2d_[i] + encode_dict['v=12']
            elif abs(res1 - res2) == 13:
                list_2d_[i] = list_2d_[i] + encode_dict['v=13']
            elif 14 <= abs(res1 - res2) < 18:
                list_2d_[i] = list_2d_[i] + encode_dict['14<=v<18']
            elif 18 <= abs(res1 - res2) < 23:
                list_2d_[i] = list_2d_[i] + encode_dict['18<=v<23']
            elif 23 <= abs(res1 - res2) < 28:
                list_2d_[i] = list_2d_[i] + encode_dict['23<=v<28']
            elif 28 <= abs(res1 - res2) < 38:
                list_2d_[i] = list_2d_[i] + encode_dict['28<=v<38']
            elif 38 <= abs(res1 - res2) < 48:
                list_2d_[i] = list_2d_[i] + encode_dict['38<=v<48']
            else:
                list_2d_[i] = list_2d_[i] + encode_dict['v>=48']
        end_time = time.time()
        print('------> metapsicov one-hot position: {time}s.'.format(time=end_time - start_time))
        return list_2d_

    def relative(
            self,
            pos,
            interval : List,
    ):
        """
        Notes
        -----
            One application scenario is to calculate the position of a residue
            relative to the boundary of a transmembrane segment (interval).
        """
        return abs(pos - min(interval)) / (max(interval) - min(interval))

    def absolute(self, pos, seq, decimal_place=10):
        return round(pos / len(seq), decimal_place)


if __name__ == "__main__":
    from pypropel.prot.sequence.Fasta import Fasta as sfasta
    from pypropel.path import to

    import tmkit as tmk
    sequence = sfasta().get(
        fasta_fpn=to("data/fasta/1aigL.fasta")
    )
    # print(sequence)

    pos_list = tmk.seq.pos_list_pair(
        len_seq=len(sequence),
        seq_sep_superior=None,
        seq_sep_inferior=0,
    )
    # print(pos_list)

    positions = tmk.seq.pos_pair(sequence=sequence, pos_list=pos_list)
    # print(positions[:3])

    window_size = 1
    win_aa_ids = tmk.seq.win_id_pair(
        sequence=sequence,
        position=positions,
        window_size=window_size,
    )
    # print(win_aa_ids)

    win_aas = tmk.seq.win_name_pair(
        sequence=sequence,
        position=positions,
        window_size=window_size,
        mids=win_aa_ids,
    )
    # print(win_aas[:3])

    features = [[] for i in range(len(win_aa_ids))]
    # print(features)
    # print(len(features))

    p = Position()
    print(positions[0][0])
    print(p.absolute(positions[0][0], sequence))

    print(p.deepconpred())

    print(p.deepconpred_(features, win_aa_ids))

    # print(p.metapsicov())

    # print(p.metapsicov_(features, win_aa_ids))