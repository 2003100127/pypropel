__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List

from pypropel.prot.feature.sequence.Composition import Composition
from pypropel.prot.feature.sequence.Length import Length


def composition(
        seq: str,
        k_spaced: int = 1,
        mode: str = 'aac',
):
    if mode == 'aac':
        return Composition(seq).aac()
    elif mode == 'dac':
        return Composition(seq).dac()
    elif mode == 'tac':
        return Composition(seq).tac()
    elif mode == 'qac':
        return Composition(seq).qac()
    elif mode == 'cksnap':
        return Composition(seq).cksnap(k=k_spaced)
    elif mode == 'aveanf':
        return Composition(seq).aveanf()
    else:
        return Composition(seq).aac()


def length(
        seq: str,
        mode: str = 'normal',
):
    if mode == 'normal':
        return Length().sequence(seq)
    elif mode == 'log':
        return Length().log(seq)
    else:
        return Length().sequence(seq)


if __name__ == "__main__":
    seq = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    print(seq)

    # print(composition(
    #     seq=seq,
    #     k_spaced=3,
    #     mode='aac',
    # ))

    # print(length(
    #     seq=seq,
    #     mode='log',
    # ))

    feature_vector = [[] for i in range(len(seq))]
    print(feature_vector)
    print(len(feature_vector))

    aac_dict = composition(
        seq=seq,
        k_spaced=3,
        mode='aac',
    )
    print(aac_dict)

    l_normal = length(
      seq=seq,
      mode='normal',
    )
    print(l_normal)

    l_log = length(
        seq=seq,
        mode='log',
    )
    print(l_log)

    for i, c in enumerate(feature_vector):
    #     print(i, c)
        feature_vector[i].append(i)
        feature_vector[i].append(seq[i])
    # print(feature_vector)

    for i, c in enumerate(feature_vector):
        feature_vector[i].append(aac_dict[feature_vector[i][1]])

    for i, c in enumerate(feature_vector):
        feature_vector[i].append(l_normal)
        feature_vector[i].append(float(l_log))

    for i, c in enumerate(feature_vector):
        if i < 32:
            feature_vector[i].append(1)
        else:
            feature_vector[i].append(0)

    # print(feature_vector)
    import numpy as np
    train_label = np.array(feature_vector)[:,-1].astype(int)
    train_feature = np.array(feature_vector)[:,2:-1].astype(float)
    # print(train_label)
    # print(train_feature)

    sd  = np.array(feature_vector)[:,2:]
    print(sd[:, np.newaxis].shape)
