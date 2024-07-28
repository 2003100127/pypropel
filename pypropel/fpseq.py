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

    # print(composition(
    #     seq=seq,
    #     k_spaced=3,
    #     mode='qac',
    # )[:3])

    print(length(
        seq=seq,
        mode='log',
    ))
