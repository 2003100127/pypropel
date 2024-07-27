__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List

from pyprocpp.prot.feature.sequence.AminoAcidProperty import AminoAcidProperty as aaprop
from pyprocpp.prot.feature.sequence.AminoAcidRepresentation import AminoAcidRepresentation as aarepr
from pyprocpp.prot.feature.Length import Length as flen


def length(
        sequence : str,
) -> int:
    return flen().sequence(sequence=sequence)


def property(
        prop_kind : str ='positive',
        prop_met : str ='Russell',
        standardize : bool =True,
) -> dict:
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
        "positive": aaprop.positive,
        "negative": aaprop.negative,
        "charged": aaprop.charged,
        "polar": aaprop.polar,
        "aliphatic": aaprop.aliphatic,
        "aromatic": aaprop.aromatic,
        "hydrophobic": aaprop.hydrophobic,
        "small": aaprop.small,
        "active": aaprop.active,
        "weight": aaprop.weight,
        "pI": aaprop.pI,
        "solubility": aaprop.solubility,
        "tm": aaprop.tm,
        "pka": aaprop.pka,
        "pkb": aaprop.pkb,
        "hydrophilicity": aaprop.hydrophilicity,
        "hydrophobicity": aaprop.hydrophobicity,
        "fet": aaprop.fet,
        "hydration": aaprop.hydration,
        "signal": aaprop.signal,
        "volume": aaprop.volume,
        "polarity": aaprop.polarity,
        "composition": aaprop.composition,
    }[prop_kind](kind=prop_met, standardize=standardize)


def onehot(
        arr_2d,
        arr_aa_names,
) -> List[List]:
    return aarepr().onehot(
        arr_2d=arr_2d,
        arr_aa_names=arr_aa_names,
    )


if __name__ == "__main__":
    from pyprocpp.prot.sequence.Fasta import Fasta as sfasta
    from pyprocpp.path import to
    import tmkit as tmk

    sequence = sfasta().get(
        fasta_fpn=to("data/fasta/1aigL.fasta")
    )
    print(sequence)

    pos_list = tmk.seq.pos_list_single(len_seq=len(sequence), seq_sep_superior=None, seq_sep_inferior=0)
    print(pos_list)

    positions = tmk.seq.pos_single(sequence=sequence, pos_list=pos_list)
    print(positions)

    win_aa_ids = tmk.seq.win_id_single(
        sequence=sequence,
        position=positions,
        window_size=1,
    )
    print(win_aa_ids)

    win_aas = tmk.seq.win_name_single(
        sequence=sequence,
        position=positions,
        window_size=1,
        mids=win_aa_ids,
    )
    print(win_aas)

    features = [[] for i in range(len(sequence))]

    print(onehot(
        arr_2d=features,
        arr_aa_names=win_aas,
    )[0])