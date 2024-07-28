__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import time
from itertools import chain
from evcouplings.couplings import CouplingsModel


class plmc:

    def __init__(self, ):
        pass

    def read(self, file_path):
        ec = CouplingsModel(file_path)
        print(ec.alphabet)
        print(ec.index_list)
        print(ec.Jij(119, 112))

    def Jij(self, list_2d, position, prot_name, file_chain, param_path):
        stime = time.time()
        pos_ = position
        ec = CouplingsModel(param_path + prot_name + file_chain + '.params')
        aa_types = len(ec.alphabet)
        print(ec.index_list)
        # print(len(ec.seq()))
        # print(ec.Jij(119, 112))
        list_2d_ = list_2d
        num_pairs = len(pos_)
        # print('Pairs number is: {}'.format(num_pairs))
        for i in range(num_pairs):
            Jij_1d = list(chain.from_iterable(ec.Jij(pos_[i][0], pos_[i][3])))
            for j in range(aa_types * aa_types):
                list_2d_[i].append(Jij_1d[j])
        # print(pos_[0])
        etime = time.time()
        print('---> plmc coupling matrix: {time}s.'.format(time=etime - stime))
        return list_2d_


if __name__ == "__main__":
    from pypropel.prot.sequence.Fasta import Fasta as sfasta
    from pypropel.path import to
    import tmkit as tmk

    sequence = sfasta().get(
        fasta_fpn=to("data/fasta/1xqfA.fasta")
    )
    print(sequence)
    print(len(sequence))

    pos_list = tmk.seq.pos_list_pair(len_seq=len(sequence), seq_sep_superior=None, seq_sep_inferior=0)
    print(pos_list)

    positions = tmk.seq.pos_pair(sequence=sequence, pos_list=pos_list)
    print(positions)

    win_aa_ids = tmk.seq.win_id_pair(
        sequence=sequence,
        position=positions,
        window_size=5,
    )
    # print(win_aa_ids)

    p = plmc()

    print(p.Jij(
        list_2d=positions,
        prot_name='1xqf',
        file_chain='A',
        position=positions,
        param_path=to('data/plmc/'),
    )[0])