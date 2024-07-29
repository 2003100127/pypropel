__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from typing import  List, Dict, Tuple

import numpy as np

from pypropel.prot.feature.alignment.MSA import MSA as msaparser
from pypropel.prot.feature.alignment.frequency.Single import Single
from pypropel.prot.feature.alignment.frequency.Pair import Pair
from pypropel.prot.feature.alignment.representation.Binary import Binary
from pypropel.prot.feature.alignment.representation.Frequency import Frequency
from pypropel.util.Console import Console


def read(
        msa_fpn : str
) -> List:
    return msaparser(
        msa_fpn=msa_fpn,
    ).read()


def split(
        msa : List
) -> List:
    msa_sp = []
    for homolog in msa:
        msa_sp.append(list(homolog))
    return msa_sp


def freq_col_sgl(
        msa: List,
):
    return Single(msa).columns()


def cnt_col_sgl(
        msa: List,
) -> Dict:
    return Single(msa).columnsByPandas()


def freq_whole_sgl(
        msa: List,
) -> Dict:
    return Single(msa).alignment()


def freq_pair(
        msa: List,
        i: int,
        j: int,
) -> Dict:
    return Pair(msa).frequency(
        i=i,
        j=j,
    )


def representation_onehot(
        msa: List,
) -> np.ndarray:
    return Binary(msa=msa).onehot()


def representation_freq(
        msa: List,
) -> np.ndarray:
    return Frequency(msa=msa).matrix()


if __name__ == "__main__":
    from pypropel.path import to

    msa = read(
        msa_fpn=to('data/msa/aln/1aijL.aln'),
    )
    # print(msa)

    print(split(msa=msa)[:1])

    # print(freq_col_sgl(msa=msa))
    # print(cnt_col_sgl(msa=msa))
    # print(freq_whole_sgl(msa=msa))
    # print(freq_pair(msa=msa, i=1, j=2))

    # oh_rep = representation_onehot(
    #     msa=msa
    # )
    # print(oh_rep)
    # print(oh_rep.shape)

    # print(representation_freq(msa=msa))