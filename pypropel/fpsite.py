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


if __name__ == "__main__":
    from pypropel.prot.sequence.Fasta import Fasta as sfasta
    from pypropel.path import to
    import tmkit as tmk

    print(property('positive'))

    sequence = sfasta().get(
        fasta_fpn=to("data/fasta/1aigL.fasta")
    )
    # print(sequence)

    pos_list = tmk.seq.pos_list_single(len_seq=len(sequence), seq_sep_superior=None, seq_sep_inferior=0)
    # print(pos_list)

    positions = tmk.seq.pos_single(sequence=sequence, pos_list=pos_list)
    # print(positions)

    win_aa_ids = tmk.seq.win_id_single(
        sequence=sequence,
        position=positions,
        window_size=1,
    )
    # print(win_aa_ids)

    win_aas = tmk.seq.win_name_single(
        sequence=sequence,
        position=positions,
        window_size=1,
        mids=win_aa_ids,
    )
    # print(win_aas)

    features = [[] for i in range(len(sequence))]
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