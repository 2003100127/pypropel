__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List, Dict

from pyprocpp.prot.feature.alignment.MSA import MSA as msaparser
from pyprocpp.prot.feature.alignment.InformationTheory import InformationTheory as itheory
from pyprocpp.prot.feature.Length import Length as flen


def read_msa(
        msa_fpn : str,
) -> List:
    return msaparser(
        msa_fpn=msa_fpn
    ).read()


def length(
        msa : str,
) -> int:
    return flen().msa(msa=msa)


def entropy(
        msa : str,
) -> Dict:
    return itheory(
        msa=msa,
    ).entropy()


def entropy_gap(
        msa : str,
        gap_thres : str = 1,
) -> Dict:
    return itheory(
        msa=msa,
    ).entropy_gap(gap_thres=gap_thres)


def mutual_information(
        msa : str,
        i : int,
        j : int,
) -> Dict:
    return itheory(
        msa=msa,
    ).mi(i=i, j=j)


if __name__ == "__main__":
    from pyprocpp.path import to

    msa = read_msa(msa_fpn=to('data/msa/1aijL.aln'))
    # print(msa)

    # print(entropy(msa=msa))
    # print(entropy_gap(msa=msa, gap_thres=100))
    print(mutual_information(msa=msa, i=1, j=2))
