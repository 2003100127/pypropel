__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Tuple

import pandas as pd
from pypropel.prot.sequence.ConvertS2M import ConvertS2M
from pypropel.prot.sequence.ConvertM2S import ConvertM2S
from pypropel.prot.feature.alignment.Convert import Convert


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


def many2single(
        fasta_fpn,
        sv_fp : str = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    cm2s = ConvertM2S(
        input_fpn=fasta_fpn,
        in_format='fasta',
        sv_fp=sv_fp,
    )
    df = cm2s.df
    df_new = cm2s.tofasta()
    return df, df_new


def msa_reformat(
        input_fpn : str,
        in_format : str,
        output_fpn : str,
        out_format : str,
):
    return Convert(
        input_fpn=input_fpn,
        in_format=in_format,
        output_fpn=output_fpn,
        out_format=out_format,
    ).reformat()


def msa2fas(
        input_fpn : str,
        output_fp : str,
):
    return Convert(
        input_fpn=input_fpn,
        output_fp=output_fp,
    ).tofasta_sgl()


if __name__ == "__main__":
    from pypropel.path import to

    # df = single2many(
    #     fasta_fp=to('data/fasta/'),
    #     prot_df=pd.DataFrame({
    #         'prot': ['1aig', '1aij', '1xqf'],
    #         'chain': ['L', 'L', 'A'],
    #     }),
    #     sv_fpn=to('data/fasta/s2m.fasta')
    # )

    # df, df_new = many2single(
    #     fasta_fpn=to('data/msa/experimental_protein.fasta'),
    #     sv_fp=to('data/msa/')
    # )

    # print(df)
    # print(df_new)

    msa_reformat(
        input_fpn=to('data/msa/a2m/ET.a2m'),
        in_format='fasta',
        output_fpn=to('data/msa/a2m/ET_converted.sto'),
        out_format='stockholm',
    )

    # msa2fas(
    #     input_fpn=to('data/msa/a2m/ET.a2m'),
    #     output_fp=to('data/msa/a2m/'),
    # )