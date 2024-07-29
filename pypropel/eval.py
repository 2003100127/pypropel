__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from typing import Tuple

import pandas as pd

from pypropel.prot.tool.isite.Dispatcher import Dispatcher


def sitewise_segment(
        prot_df : pd.DataFrame,
        dist_fp : str,
        dist_limit : float,
        tool_fp : str,
        tool : str,
        pdb_fp : str,
        topo_fp : str,
        xml_fp : str,
        fasta_fp : str,
        segment : str,
        sort : int,
        sv_fp=None,
):
    return Dispatcher().segment(
        prot_df=prot_df,
        dist_fp=dist_fp,
        dist_limit=dist_limit,
        tool_fp=tool_fp,
        tool=tool,
        pdb_fp=pdb_fp,
        topo_fp=topo_fp,
        xml_fp=xml_fp,
        fasta_fp=fasta_fp,
        segment=segment,
        sort=sort,
        sv_fp=sv_fp,
    )


if __name__ == "__main__":
    from pypropel.path import to

    prot_df = pd.DataFrame({
        'prot': ['1aij', ],
        'chain': ['L', ],
    })

    sitewise_segment(
        prot_df=prot_df,
        dist_fp=to('data/pdb/complex/pdbtm/'),
        dist_limit=6.,
        tool_fp=to('data/isite/deeptminter/'),
        tool='tma300',
        pdb_fp=to('data/pdb/pdbtm/'),
        topo_fp=to('data/phobius/'),
        xml_fp=to('data/xml/'),
        fasta_fp=to('data/fasta/'),
        segment='pdbtm_tmh',
        sort=1,
        sv_fp=to('data/'),
    )