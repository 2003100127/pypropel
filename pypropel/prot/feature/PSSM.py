__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import time
import numpy as np
from pypropel.util.Console import Console


class PSSM:

    def __init__(
            self,
            verbose: bool = True,
    ):
        self.console = Console()
        self.console.verbose = verbose

    def blast(
            self,
            blast_fpn,
    ):
        """"""
        from pypropel.prot.sequence.Symbol import Symbol
        aa_universal = Symbol().single(gap=False, universal=True)
        aa_vocabulary = {i: 0 for i in Symbol().single(gap=False, universal=False)}
        pssm_ = {}
        f = open(blast_fpn, 'r')
        c = 0
        for line in f.readlines():
            split_line = line.replace('-', ' -').split()
            if (len(split_line) in (44, 22)) and (split_line[0] not in ('#', 'Last')):
                c = c + 1
                pssm_universal = {
                    aa_universal[i]: -float(e) for i, e in enumerate(split_line[2:22])
                }
                for k, v in pssm_universal.items():
                    aa_vocabulary[k] = v
                # print(pssm_universal)
                # print(aa_vocabulary)
                pssm_[c] = [v for v in aa_vocabulary.values()]
                # print(pssm_)
        f.close()
        return pssm_

    def hhm(
            self,
            hhm_fpn,
    ):
        f = open(hhm_fpn)
        line = f.readline()
        while line[0] != '#':
            line = f.readline()
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        line = f.readline()
        c = 0
        hhm_ = {}
        while line[0: 2] != '//':
            c += 1
            lineinfo = line.split()
            probs_ = [2 ** (-float(lineinfo[i]) / 1000) if lineinfo[i] != '*' else 0. for i in range(2, 22)]
            probs_ = np.array(probs_)[np.newaxis, :]
            # print(probs_)
            line = f.readline()
            lineinfo = line.split()
            extras_ = [2 ** (-float(lineinfo[i]) / 1000) if lineinfo[i] != '*' else 0. for i in range(0, 10)]
            extras_ = np.array(extras_)[np.newaxis, :]
            # print(extras_)
            combine = np.squeeze(np.concatenate((probs_, extras_), axis=1), axis=0)
            # print(combine)
            line = f.readline()
            assert len(line.strip()) == 0
            line = f.readline()
            hhm_[c] = combine
        return hhm_

    def blast_(
            self,
            pssm,
            list_2d,
            window_aa_ids,
    ):
        list_2d_ = list_2d
        for i, aa_win_ids in enumerate(window_aa_ids):
            for j in aa_win_ids:
                # print(j)
                for k in j:
                    # print(k)
                    if k is None:
                        list_2d_[i] = list_2d_[i] + [0 for i in range(20)]
                    else:
                        list_2d_[i] = list_2d_[i] + pssm[k]
        return list_2d_

    def hhm_(
            self,
            hhm,
            list_2d,
            window_aa_ids,
    ):
        list_2d_ = list_2d
        for i, aa_win_ids in enumerate(window_aa_ids):
            for j in aa_win_ids:
                # print(j)
                for k in j:
                    # print(k)
                    if k is None:
                        list_2d_[i] = list_2d_[i] + [0 for i in range(30)]
                    else:
                        list_2d_[i] = list_2d_[i] + hhm[k].tolist()
        return list_2d_


if __name__ == "__main__":
    from pypropel.prot.sequence.Fasta import fasta as sfasta
    from pypropel.prot.feature.window.Pair import pair as pwindow
    from pypropel.prot.position.Fasta import fasta as pfasta
    from pypropel.prot.position.scenario.Length import length as lscenario

    INIT = {
        'prot_name': '1aig',
        'file_chain': 'L',
        'msa_path': to('data/protein/msa/tm_alpha_n165/'),
        'fasta_path': to('data/protein/fasta/tm_alpha_n165/'),
        'pssm_path': to('data/protein/fasta/tm_alpha_n165/'),
        'hhm_path': to('data/protein/fasta/tm_alpha_n165/'),
    }

    window_size = 0
    sequence = sfasta().get(
        fasta_path=INIT['fasta_path'],
        fasta_name=INIT['prot_name'],
        file_chain=INIT['file_chain']
    )

    # ### ++++++++++++++++++ start ppi ++++++++++++++++++++
    # length_pos_list = lscenario().toSingle(len(sequence))
    # positions = pfasta(sequence).single(length_pos_list)
    # window_aa_ids = swindow(
    #     sequence=sequence,
    #     position=positions,
    #     window_size=window_size,
    # ).aaid()
    # print(window_aa_ids)
    # ### +++++++++++++++++++ end ppi +++++++++++++++++++

    ### ++++++++++++++++++ start rrc ++++++++++++++++++++
    length_pos_list = lscenario().toPair(len(sequence))
    positions = pfasta(sequence).pair(length_pos_list)
    window_aa_ids = pwindow(
        sequence=sequence,
        position=positions,
        window_size=window_size,
    ).aaid()
    print(window_aa_ids)
    ### +++++++++++++++++++ end rrc +++++++++++++++++++

    # features_1d_in = [[] for i in range(len(sequence))]
    features_2d_in = length_pos_list

    # msa = msaparser(INIT['msa_path'] + INIT['prot_name'] + INIT['file_chain'] + '.aln').read()

    # p = evolutionaryProfile(msa)

    # print(p.scheme1())

    # print(p.scheme2())

    # print(p.scheme1_(features, window_aa_ids))

    # print(p.scheme2_(features, window_aa_ids))

    p = PSSM()

    # print(p.pssm(
    #     pssm_fpn=to('data/pssm/1aigL.pssm'),
    # ))

    print(p.hhm(
        hhm_fpn=to('data/hhm/1aigL.hhm'),
    ))

    print(p.blast__(
        pssm_path=to('data/predictor/profile/pssm/tm_alpha_n165/'),
        prot_name=INIT['prot_name'],
        file_chain=INIT['file_chain'],
        list_2d=features_2d_in,
        window_aa_ids=window_aa_ids,
    ))

    print(p.hhm__(
        hhm_path=to('data/predictor/profile/hhm/tm_alpha_n165/'),
        prot_name=INIT['prot_name'],
        file_chain=INIT['file_chain'],
        list_2d=features_2d_in,
        window_aa_ids=window_aa_ids,
    ))