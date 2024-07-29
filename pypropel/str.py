__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from typing import  List, Dict

import pandas as pd

from pypropel.prot.sequence.PDB import PDB
from pypropel.prot.structure.convert.ToFasta import ToFasta
from pypropel.prot.structure.chain.Format import Format
from pypropel.prot.structure.chain.Splitter import Splitter
from pypropel.prot.structure.hetatm.Remove import Remove
from pypropel.prot.structure.distance.isite.heavy.AllAgainstAll import AllAgainstAll
from pypropel.util.Console import Console


def read(
        pdb_path,
        pdb_name,
        file_chain,
        seq_chain,
):
    return PDB(
        pdb_path=pdb_path,
        pdb_name=pdb_name,
        file_chain=file_chain,
        seq_chain=seq_chain,
    ).chain()


def chains(
        pdb_fp,
        pdb_name,
):
    return AllAgainstAll(
        pdb_fp=pdb_fp,
        pdb_name=pdb_name,
    ).chains()


def tofasta(
        prot_df: pd.DataFrame,
        sv_fp: str,
        pdb_path: str,
):
    return ToFasta(
        prot_df=prot_df,
        sv_fp=sv_fp,
    ).frompdb(
        pdb_path=pdb_path,
    )


def del_end(
        prot_df: pd.DataFrame,
        sv_fp: str,
        pdb_path: str,
):
    return Format(
        prot_df=prot_df,
        sv_fp=sv_fp,
    ).del_END_frompdb(
        pdb_path=pdb_path,
    )


def split_cplx_to_sgl(
        prot_df: pd.DataFrame,
        pdb_path: str,
        sv_fp: str,
):
    return Splitter(
        prot_df=prot_df,
        pdb_path=pdb_path,
        sv_fp=sv_fp,
    ).pdb_per_chain()


def remove_hetatm(
        prot_df: pd.DataFrame,
        pdb_path: str,
        sv_fp: str,
):
    return Remove(
        prot_df=prot_df,
    ).biopython(
        pdb_path=pdb_path,
        sv_fp=sv_fp,
    )


if __name__ == "__main__":
    from pypropel.path import to

    # print(read(
    #     pdb_path=to('data/pdb/pdbtm/'),
    #     pdb_name='1aij',
    #     file_chain='L',
    #     seq_chain='L',
    # ))

    print(chains(
        pdb_fp=to('data/pdb/complex/pdbtm/'),
        pdb_name='1aij',
    ))

    prot_df = pd.DataFrame({
        'prot': ['1aig', '1aij', '1xqf'],
        'chain': ['L', 'L', 'A'],
    })

    # print(tofasta(
    #     prot_df,
    #     sv_fp=to('data/'),
    #     pdb_path=to('data/pdb/pdbtm/'),
    # ))

    # print(del_end(
    #     prot_df,
    #     sv_fp=to('data/'),
    #     pdb_path=to('data/pdb/pdbtm/'),
    # ))

    # print(split_cplx_to_sgl(
    #     prot_df=prot_df,
    #     pdb_path=to('data/pdb/complex/pdbtm/'),
    #     sv_fp=to('data/'),
    # ))

    # print(remove_hetatm(
    #     prot_df=prot_df,
    #     pdb_path=to('data/pdb/complex/pdbtm/'),
    #     sv_fp=to('data/pdb/'),
    # ))