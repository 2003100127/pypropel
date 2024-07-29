__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from typing import  List, Dict

import pandas as pd

from pypropel.prot.sequence.Fasta import Fasta
from pypropel.prot.sequence.IsEmpty import IsEmpty
from pypropel.prot.sequence.IsMatch import IsMatch
from pypropel.util.Console import Console


def read(
        fasta_fpn: str
):
    return Fasta().get(fasta_fpn)


def save(
        list_2d: list[list[str]],
        sv_fp: str,
):
    return Fasta().save(
        list_2d=list_2d,
        sv_fp=sv_fp
    )


def is_empty(
    prot_df: pd.DataFrame,
    fasta_fp: str,
    sv_empty_fp: str,
):
    return IsEmpty(
        prot_df=prot_df,
        fasta_fp=fasta_fp,
        sv_empty_fp=sv_empty_fp,
    ).fasta()


def is_match(
        prot_df : pd.DataFrame,
        kind : str,
        sv_mismatch_fp : str,
        fasta_path: str = None,
        pdb_path: str = None,
        xml_path: str = None,
):

    return IsMatch(
        prot_df=prot_df,
        fasta_path=fasta_path,
        pdb_path=pdb_path,
        xml_path=xml_path,
        kind=kind,
        sv_mismatch_fp=sv_mismatch_fp,
    ).execute()


if __name__ == "__main__":
    from pypropel.path import to


    seq = read(
        fasta_fpn=to('data/fasta/1aigL.fasta'),
    )
    print(seq)

    # print(save(
    #     list_2d=[
    #         ['1aigL-new', seq],
    #     ],
    #     sv_fp=to('data/fasta/'),
    # ))

    # prot_df = pd.DataFrame({
    #     'prot': ['1aig', '1aij', '1xqf'],
    #     'chain': ['L', 'L', 'A'],
    # })

    # print(is_empty(
    #     prot_df=prot_df,
    #     fasta_fp=to('data/fasta/'),
    #     sv_empty_fp=to('data/'),
    # ))

    # print(is_match(
    #     prot_df=prot_df,
    #     fasta_path=to('data/fasta/'),
    #     # pdb_path=to('data/pdb/pdbtm/'),
    #     xml_path=to('data/xml/'),
    #     kind='fasta<->xml',
    #     sv_mismatch_fp=to('data/'),
    # ))