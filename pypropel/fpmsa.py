__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List, Dict

from pypropel.prot.feature.alignment.Composition import Composition as compo
from pypropel.prot.feature.PSSM import PSSM
from pypropel.prot.feature.alignment.InformationTheory import InformationTheory as itheory
from pypropel.prot.feature.alignment.Conservation import Conservation
from pypropel.prot.feature.sequence.Length import Length as flen


def length(
        msa : List,
) -> int:
    return flen().msa(msa=msa)


def entropy(
        msa : List,
) -> Dict:
    return itheory(
        msa=msa,
    ).entropy()


def entropy_gap(
        msa : List,
        gap_thres : int = 1,
) -> Dict:
    return itheory(
        msa=msa,
    ).entropy_gap(gap_thres=gap_thres)


def conservation_custom(
        ent_dict: Dict,
) -> Dict:
    return Conservation().get(
        ent_dict
    )


def mutual_information(
        msa : List,
        i : int,
        j : int,
) -> Dict:
    return itheory(
        msa=msa,
    ).mi(i=i, j=j)


def composition(
        msa : List,
        mode: str = 'aac'
):
    if mode == 'aac':
        return compo(msa=msa).aac()
    else:
        return compo(msa=msa).aac()


def pssm(
        fpn : str = None,
        msa: List = None,
        mode : str = 'blast'
):
    if mode == 'blast':
        return PSSM().blast(
            blast_fpn=fpn
        )
    elif mode == 'hhm':
        return PSSM().hhm(
            hhm_fpn=fpn
        )
    elif mode == 'ep':
        return compo(msa=msa).evolutionary_profile()
    elif mode == 'ep_norm':
        return compo(msa=msa).evolutionary_profile_norm()
    else:
        return PSSM().blast(
            blast_fpn=fpn
        )


def jsd(
        fpn : str = None,
        mode : str = 'standalone'
):
    if mode == 'standalone':
        return Conservation().jsd(
            jsd_fpn=fpn
        )


def consurf(
        fpn : str = None,
        mode : str = 'v1'
):
    if mode == 'v1':
        return Conservation().consurf_v1(
            consurf_fpn=fpn
        )


if __name__ == "__main__":
    from pypropel.prot.feature.alignment.MSA import MSA as msaparser
    from pypropel.path import to

    msa = msaparser(msa_fpn=to('data/msa/aln/1aijL.aln')).read()
    # print(msa)

    # print(length(msa=msa))

    ent_dict = entropy(msa=msa)
    # print(ent_dict)

    # ent_dict = entropy_gap(msa=msa, gap_thres=100)
    # print(ent_dict)

    # print(conservation_custom(
    #     ent_dict=ent_dict,
    # ))

    # print(mutual_information(msa=msa, i=1, j=2))

    print(composition(
        msa=msa,
        mode='aac',
    ))

    # print(pssm(
    #     # fpn=to('data/pssm/1aigL.pssm'),
    #     # mode='blast',
    #     fpn=to('data/hhm/1aigL.hhm'),
    #     mode='hhm',
    # ))

    # print(jsd(
    #     fpn=to('data/conservation/jsd/SR24_CtoU/CLEC2B_LOC113845378.jsd'),
    #     mode='standalone',
    # ))

    # print(consurf(
    #     fpn=to('data/conservation/consurf/E.consurf'),
    #     mode='v1'
    # ))