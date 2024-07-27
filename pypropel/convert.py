__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List, Dict

import pandas as pd
from pyprocpp.prot.sequence.ConvertS2M import ConvertS2M
from pyprocpp.prot.sequence.ConvertM2S import ConvertM2S
from pyprocpp.prot.feature.alignment.Convert import Convert


def single2many(
        fasta_fp,
        prot_df,
        sv_fpn : str = None
) -> pd.DataFrame:
    cs2m = ConvertS2M(
        fasta_fp=fasta_fp,
        prot_df=prot_df,
    )
    seqs = cs2m.integrate_seq()
    if sv_fpn:
        cs2m.save(
            seqs,
            sv_fpn=sv_fpn
        )
    return pd.DataFrame(seqs)


if __name__ == "__main__":
    from pyprocpp.path import to

    df = single2many(
        fasta_fp=to('data/fasta/'),
        prot_df=pd.DataFrame({
            'prot': ['1aig', '1aij', '1xqf'],
            'chain': ['L', 'L', 'A'],
        }),
        sv_fpn=to('data/fasta/s2m.fasta')
    )