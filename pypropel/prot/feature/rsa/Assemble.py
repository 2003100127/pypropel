__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import time
import numpy as np
from pypropel.util.Console import Console


class Assemble:

    def __init__(
            self,
            verbose: bool = True,

    ):
        self.console = Console()
        self.console.verbose = verbose

    def accpro(
            self,
            df_accpro,
            list_2d,
            window_aa_ids,
            mode='single',
    ):
        start_time = time.time()
        accpro = df_accpro.values.tolist()
        print(df_accpro)
        window_aa_ids_ = window_aa_ids if mode == 'single' else [i[0] + i[1] for i in window_aa_ids]
        list_2d_ = list_2d
        for i, aa_win_ids in enumerate(window_aa_ids_):
            # print(i)
            # print(aa_win_ids)
            for j in aa_win_ids:
                # print(aa_win_ids)
                if j is None:
                    list_2d_[i] = list_2d_[i] + np.zeros(2).tolist()
                else:
                    if accpro[j - 1][0] == 'e':
                        list_2d_[i] = list_2d_[i] + [0, 1]
                    else:
                        list_2d_[i] = list_2d_[i] + [1, 0]
        end_time = time.time()
        self.console.print('=========>ACCpro solvent: {time}s.'.format(time=end_time - start_time))
        return list_2d_

    def accpro20(
            self,
            df_accpro20,
            list_2d,
            window_aa_ids,
            mode='single',
    ):
        start_time = time.time()
        accpro20 = df_accpro20.values.tolist()
        # print(accpro20)
        window_aa_ids_ = window_aa_ids if mode == 'single' else [i[0] + i[1] for i in window_aa_ids]
        list_2d_ = list_2d
        for i, aa_win_ids in enumerate(window_aa_ids_):
            # print(i)
            for j in aa_win_ids:
                # print(aa_win_ids)
                if j is None:
                    list_2d_[i] = list_2d_[i] + np.zeros(20).tolist()
                else:
                    if accpro20[j - 1][0] == 0.0:
                        list_2d_[i] = list_2d_[i] + [1 if i == 19 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.05:
                        list_2d_[i] = list_2d_[i] + [1 if i == 18 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.1:
                        list_2d_[i] = list_2d_[i] + [1 if i == 17 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.15:
                        list_2d_[i] = list_2d_[i] + [1 if i == 16 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.2:
                        list_2d_[i] = list_2d_[i] + [1 if i == 15 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.25:
                        list_2d_[i] = list_2d_[i] + [1 if i == 14 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.3:
                        list_2d_[i] = list_2d_[i] + [1 if i == 13 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.35:
                        list_2d_[i] = list_2d_[i] + [1 if i == 12 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.4:
                        list_2d_[i] = list_2d_[i] + [1 if i == 11 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.45:
                        list_2d_[i] = list_2d_[i] + [1 if i == 10 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.5:
                        list_2d_[i] = list_2d_[i] + [1 if i == 9 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.55:
                        list_2d_[i] = list_2d_[i] + [1 if i == 8 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.6:
                        list_2d_[i] = list_2d_[i] + [1 if i == 7 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.65:
                        list_2d_[i] = list_2d_[i] + [1 if i == 6 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.7:
                        list_2d_[i] = list_2d_[i] + [1 if i == 5 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.75:
                        list_2d_[i] = list_2d_[i] + [1 if i == 4 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.20:
                        list_2d_[i] = list_2d_[i] + [1 if i == 3 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.205:
                        list_2d_[i] = list_2d_[i] + [1 if i == 2 else 0 for i in range(20)]
                    elif accpro20[j - 1][0] == 0.9:
                        list_2d_[i] = list_2d_[i] + [1 if i == 1 else 0 for i in range(20)]
                    else:
                        list_2d_[i] = list_2d_[i] + [1 if i == 0 else 0 for i in range(20)]
        end_time = time.time()
        self.console.print('=========>ACCpro20 solvent: {time}s.'.format(time=end_time - start_time))
        return list_2d_

    def solvpred(
            self,
            df_solvpred,
            list_2d,
            window_aa_ids,
            mode='single',
    ):
        start_time = time.time()
        solvpred = df_solvpred.values.tolist()
        # print(solvpred)
        window_aa_ids_ = window_aa_ids if mode == 'single' else [i[0] + i[1] for i in window_aa_ids]
        list_2d_ = list_2d
        for i, aa_win_ids in enumerate(window_aa_ids_):
            # print(i)
            for j in aa_win_ids:
                # print(aa_win_ids)
                if j is None:
                    list_2d_[i].append(0)
                else:
                    list_2d_[i].append(solvpred[j-1][2])
        end_time = time.time()
        print('=========>solvpred solvent: {time}s.'.format(time=end_time - start_time))
        return list_2d_


if __name__ == "__main__":
    from pypropel.prot.sequence.Fasta import Fasta as sfasta
    from pypropel.path import to
    import tmkit as tmk

    sequence = sfasta().get(
        fasta_fpn=to("data/fasta/1aigL.fasta")
    )
    print(sequence)

    pos_list = tmk.seq.pos_list_pair(len_seq=len(sequence), seq_sep_superior=None, seq_sep_inferior=0)
    # print(pos_list)

    positions = tmk.seq.pos_pair(sequence=sequence, pos_list=pos_list)
    # print(positions)

    window_size = 0
    win_aa_ids = tmk.seq.win_id_pair(
        sequence=sequence,
        position=positions,
        window_size=window_size,
    )
    print(win_aa_ids)

    features_1d_in = [[] for i in range(len(sequence))]
    features_2d_in = positions

    from pypropel.prot.feature.rsa.Reader import Reader as a11yreader

    p = Assemble()
    # df_accpro = a11yreader().accpro(
    #     accpro_path=to('data/accessibility/accpro/'),
    #     prot_name='1aig',
    #     file_chain='L'
    # )
    # print(p.accpro(
    #     df_accpro=df_accpro,
    #     list_2d=positions,
    #     window_aa_ids=win_aa_ids,
    #     mode='pair'
    # ))

    # df_accpro20 = a11yreader().accpro20(
    #     accpro20_path=to('data/accessibility/accpro20/'),
    #     prot_name='1aig',
    #     file_chain='L'
    # )
    # print(p.accpro20(
    #     df_accpro20,
    #     list_2d=positions,
    #     window_aa_ids=win_aa_ids,
    #     mode='pair'
    # ))

    df_solvpred = a11yreader().solvpred(
        solvpred_fp=to('data/accessibility/solvpred/'),
        prot_name='1aig',
        file_chain='L',
    )
    print(p.solvpred(
        df_solvpred=df_solvpred,
        list_2d=positions,
        window_aa_ids=win_aa_ids,
        mode='pair'
    ))