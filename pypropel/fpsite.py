__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List, Dict

from pypropel.prot.feature.sequence.AminoAcidProperty import AminoAcidProperty as aaprop
from pypropel.prot.feature.sequence.AminoAcidRepresentation import AminoAcidRepresentation as aarepr
from pypropel.prot.feature.sequence.Position import Position
from pypropel.prot.feature.rsa.Reader import Reader as rsareader
from pypropel.prot.feature.ss.Reader import Reader as ssreader


def property(
        prop_kind : str ='positive',
        prop_met : str ='Russell',
        standardize : bool =True,
) -> Dict:
    """
    An amino acid's property

    Parameters
    ----------
    prop_kind
        an amino acid's property kind
    prop_met
        method from which a property is derived,
    standalize
        if standardization

    Returns
    -------

    """
    return {
        "positive": aaprop().positive,
        "negative": aaprop().negative,
        "charged": aaprop().charged,
        "polar": aaprop().polar,
        "aliphatic": aaprop().aliphatic,
        "aromatic": aaprop().aromatic,
        "hydrophobic": aaprop().hydrophobic,
        "small": aaprop().small,
        "active": aaprop().active,
        "weight": aaprop().weight,
        "pI": aaprop().pI,
        "solubility": aaprop().solubility,
        "tm": aaprop().tm,
        "pka": aaprop().pka,
        "pkb": aaprop().pkb,
        "hydrophilicity": aaprop().hydrophilicity,
        "hydrophobicity": aaprop().hydrophobicity,
        "fet": aaprop().fet,
        "hydration": aaprop().hydration,
        "signal": aaprop().signal,
        "volume": aaprop().volume,
        "polarity": aaprop().polarity,
        "composition": aaprop().composition,
    }[prop_kind](kind=prop_met, standardize=standardize)


def onehot(
        arr_2d,
        arr_aa_names,
) -> List[List]:
    return aarepr().onehot(
        arr_2d=arr_2d,
        arr_aa_names=arr_aa_names,
    )


def pos_abs_val(
        pos : int,
        seq : str,
):
    return Position().absolute(
        pos=pos,
        seq=seq,
    )

def pos_rel_val(
        pos : int,
        interval : List,
):
    return Position().relative(
        pos=pos,
        interval=interval,
    )


def deepconpred():
    return Position().deepconpred()


def metapsicov():
    return Position().metapsicov()


def rsa_solvpred(
        solvpred_fp,
        prot_name,
        file_chain,
):
    return rsareader().solvpred(
        solvpred_fp=solvpred_fp,
        prot_name=prot_name,
        file_chain=file_chain,
    )


def rsa_accpro(
        accpro_fp,
        prot_name,
        file_chain,
):
    return rsareader().accpro(
        accpro_fp=accpro_fp,
        prot_name=prot_name,
        file_chain=file_chain,
    )


def rsa_accpro20(
        accpro20_fp,
        prot_name,
        file_chain,
):
    return rsareader().accpro20(
        accpro20_fp=accpro20_fp,
        prot_name=prot_name,
        file_chain=file_chain,
    )


def ss_spider3(
        spider3_path,
        prot_name,
        file_chain,
):
    return ssreader().spider3(
        spider3_path=spider3_path,
        prot_name=prot_name,
        file_chain=file_chain,
    )


def ss_spider3_ss(
        spider3_path,
        prot_name,
        file_chain,
        sv_fp,
):
    return ssreader().spider3_to_ss(
        spider3_path=spider3_path,
        prot_name=prot_name,
        file_chain=file_chain,
        sv_fp=sv_fp,
    )


def ss_psipred(
        psipred_path,
        prot_name,
        file_chain,
        kind='ss'
):
    """

    Parameters
    ----------
    psipred_path
    prot_name
    file_chain
    kind
        1. ss;
        2. ss2;
        3. horiz

    Returns
    -------

    """
    if kind == 'ss':
        return ssreader().psipred(
            psipred_ss_path=psipred_path,
            prot_name=prot_name,
            file_chain=file_chain,
        )
    if kind == 'ss2':
        return ssreader().psipred(
            psipred_ss2_path=psipred_path,
            prot_name=prot_name,
            file_chain=file_chain,
        )
    if kind == 'horiz':
        return ssreader().psipred(
            psipred_horiz_path=psipred_path,
            prot_name=prot_name,
            file_chain=file_chain,
        )
    else:
        return ssreader().psipred(
            psipred_ss_path=psipred_path,
            prot_name=prot_name,
            file_chain=file_chain,
        )


def ss_sspro(
        sspro_path,
        prot_name,
        file_chain,
):
    return ssreader().sspro(
        sspro_path=sspro_path,
        prot_name=prot_name,
        file_chain=file_chain,
    )


def ss_sspro8(
        sspro8_path,
        prot_name,
        file_chain,
):
    return ssreader().sspro8(
        sspro8_path=sspro8_path,
        prot_name=prot_name,
        file_chain=file_chain,
    )


if __name__ == "__main__":
    from pypropel.prot.sequence.Fasta import Fasta as sfasta
    from pypropel.path import to
    # import tmkit as tmk

    # print(property('positive'))
    #
    # sequence = sfasta().get(
    #     fasta_fpn=to("data/fasta/1aigL.fasta")
    # )
    # # print(sequence)
    #
    # pos_list = tmk.seq.pos_list_single(len_seq=len(sequence), seq_sep_superior=None, seq_sep_inferior=0)
    # # print(pos_list)
    #
    # positions = tmk.seq.pos_single(sequence=sequence, pos_list=pos_list)
    # # print(positions)
    #
    # win_aa_ids = tmk.seq.win_id_single(
    #     sequence=sequence,
    #     position=positions,
    #     window_size=1,
    # )
    # # print(win_aa_ids)
    #
    # win_aas = tmk.seq.win_name_single(
    #     sequence=sequence,
    #     position=positions,
    #     window_size=1,
    #     mids=win_aa_ids,
    # )
    # # print(win_aas)
    #
    # features = [[] for i in range(len(sequence))]
    # print(features)
    # print(len(features))

    # print(onehot(
    #     arr_2d=features,
    #     arr_aa_names=win_aas,
    # )[0])

    # print(pos_abs_val(
    #     pos=positions[0][0],
    #     seq=sequence,
    # ))
    #
    # print(pos_rel_val(
    #     pos=positions[0][0],
    #     interval=[4, 10],
    # ))
    #
    # print(deepconpred())
    #
    # print(metapsicov())

    # print(rsa_solvpred(
    #     solvpred_fp=to('data/accessibility/solvpred/'),
    #     prot_name='1aig',
    #     file_chain='L',
    # ))

    # print(rsa_accpro(
    #     accpro_fp=to('data/accessibility/accpro/'),
    #     prot_name='1aig',
    #     file_chain='L',
    # ))

    # print(rsa_accpro20(
    #     accpro20_fp=to('data/accessibility/accpro20/'),
    #     prot_name='1aig',
    #     file_chain='L',
    # ))

    # print(ss_psipred(
    #     psipred_path=to('data/ss/psipred/'),
    #     prot_name='1aig',
    #     file_chain='L',
    #     kind='ss', # horiz, ss, ss2
    # ))

    # print(ss_sspro(
    #     sspro_path=to('data/ss/sspro/'),
    #     prot_name='1aig',
    #     file_chain='L'
    # ))

    # print(ss_sspro8(
    #     sspro8_path=to('data/ss/sspro8/'),
    #     prot_name='1aig',
    #     file_chain='L'
    # ))

    # print(ss_spider3(
    #     spider3_path=to('data/ss/spider3/'),
    #     prot_name='E',
    #     file_chain=''
    # ))

    # print(ss_spider3_ss(
    #     spider3_path=to('data/ss/spider3/'),
    #     prot_name='E',
    #     file_chain='',
    #     sv_fp=to('data/ss/spider3/'),
    # ))


    seq = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"
    print(seq)

    feature_vector = [[] for i in range(len(seq))]
    print(feature_vector)
    print(len(feature_vector))

    print(property(
        prop_met='Hopp',
        prop_kind='hydrophilicity'
    ))

