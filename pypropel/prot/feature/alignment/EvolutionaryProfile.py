__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np
import time
from pyprocpp.prot.feature.alignment.frequency.Single import Single as fs
from protein.tool.profile.Blast import blast
from protein.tool.profile.HHblits import hhblits


class evolutionaryProfile(object):
    """
    .. @methods:
       ---------
       getEntireEvolutionaryProfile().

    .. @summary:
       ---------
       multiSequenceRetrievalTools class mainly use tools of multiple sequence alignment.

    .. @since:
       -------
       It was introduced since v1.0.
    """

    def __init__(self, msa=None):
        """
        .. @summary:
           ---------
           parameters initialization.

        :param msa_path:
        :param msa_name:
        """
        if msa is not None:
            self.msa = msa
            self.fs = fs(self.msa)

    def scheme1(self):
        """"""
        freq_column, _ = self.fs.columns()
        # print(freq_column)
        fre_msa = self.fs.alignment()
        fraction = {}
        for aa in self.freq.aa_alphabet:
            if fre_msa[aa] == 0:
                fraction[aa] = freq_column[aa] * 0
            else:
                fraction[aa] = freq_column[aa] / fre_msa[aa]
        return fraction

    def scheme1_(self, list_2d, window_aa_ids):
        start_time = time.time()
        list_2d_ = list_2d
        ep = self.scheme1()
        # print(ep)
        for i, aa_win_ids in enumerate(window_aa_ids):
            for j in aa_win_ids:
                for aa in self.freq.aa_alphabet:
                    if j is None:
                        list_2d_[i].append(0)
                    else:
                        list_2d_[i].append(ep[aa][j-1])
        # print(len(list_2d_[0]))
        # for i in list_2d_:
        #     if len(i) != 189:
        #         print(len(i))
        end_time = time.time()
        print('------> evolutionary profile: {time}s.'.format(time=end_time - start_time))
        return list_2d_

    def scheme2(self):
        """
        get by Hoenigschmid
        :return:
        """
        fraction = self.scheme1()
        for aa in self.freq.aa_alphabet:
            for j in range(len(fraction['A'])):
                if fraction[aa][j] == 0:
                    fraction[aa][j] = 1
                else:
                    part = -np.log(fraction[aa][j])
                    fraction[aa][j] = 1 / (1 + np.exp(part))
        return fraction

    def scheme2_(self, list_2d, window_aa_ids):
        start_time = time.time()
        list_2d_ = list_2d
        ep = self.scheme2()
        for i, aa_win_ids in enumerate(window_aa_ids):
            for j in aa_win_ids:
                for aa in self.freq.aa_alphabet:
                    if j is None:
                        list_2d_[i].append(0)
                    else:
                        list_2d_[i].append(ep[aa][j - 1])
        # print(len(list_2d_[0]))
        # for i in list_2d_:
        #     if len(i) != 189:
        #         print(len(i))
        end_time = time.time()
        print('------> evolutionary profile: {time}s.'.format(time=end_time - start_time))
        return list_2d_

    def blast__(self, pssm_path, prot_name, file_chain, list_2d, window_aa_ids):
        profile = blast().pssm(
            pssm_path=pssm_path,
            prot_name=prot_name,
            file_chain=file_chain,
        )
        # print(profile)
        list_2d_ = list_2d
        for i, aa_win_ids in enumerate(window_aa_ids):
            for j in aa_win_ids:
                # print(j)
                for k in j:
                    # print(k)
                    if k is None:
                        list_2d_[i] = list_2d_[i] + [0 for i in range(20)]
                    else:
                        list_2d_[i] = list_2d_[i] + profile[k]
        return list_2d_

    def hhm__(self, hhm_path, prot_name, file_chain, list_2d, window_aa_ids):
        profile = hhblits().hhm(
            hhm_path=hhm_path,
            prot_name=prot_name,
            file_chain=file_chain,
        )
        # print(profile)
        list_2d_ = list_2d
        for i, aa_win_ids in enumerate(window_aa_ids):
            for j in aa_win_ids:
                # print(j)
                for k in j:
                    # print(k)
                    if k is None:
                        list_2d_[i] = list_2d_[i] + [0 for i in range(30)]
                    else:
                        list_2d_[i] = list_2d_[i] + profile[k].tolist()
        return list_2d_


if __name__ == "__main__":
    from protein.feature.window.Pair import pair as pwindow
    from protein.sequence.Fasta import fasta as sfasta
    from protein.position.Fasta import fasta as pfasta
    from protein.position.scenario.Length import length as lscenario

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

    p = evolutionaryProfile()

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