__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from typing import Dict

import pandas as pd

from pypropel.prot.feature.structure.DSSP import DSSP
from pypropel.prot.feature.structure.Threedi import Threedi


def threedi(
        prot_name,
        prot_chain,
        pdb_fp,
        mode='chain',
) -> Dict:
    return Threedi(
        prot_name=prot_name,
        prot_chain=prot_chain,
        pdb_fp=pdb_fp,
    ).encode(mode=mode)


def dssp_rsa_run(
        prot_name,
        prot_chain,
        pdb_fp,
        sv_fp,
) -> pd.DataFrame:
    return DSSP(
        prot_name=prot_name,
        prot_chain=prot_chain,
    ).run(
        pdb_fp=pdb_fp,
        sv_fp=sv_fp,
    )


def dssp_rsa_access(
        prot_name,
        prot_chain,
        rsa_fp,
) -> pd.DataFrame:
    return DSSP(
        prot_name=prot_name,
        prot_chain=prot_chain,
    ).access(
        rsa_fp=rsa_fp,
    )


if __name__ == "__main__":
    from pypropel.path import to
    import pandas as pd

    # prot_df = pd.DataFrame({
    #     'prot': ['1aij', '1aig', '1xqf', ],
    #     'chain': ['L', 'L', 'A', ],
    # })
    # for i in prot_df.index:
    #     print('No.{}: protein: {} chain: {}'.format(i + 1, prot_df.loc[i, 'prot'], prot_df.loc[i, 'chain']))
    #     threedi_dict = threedi(
    #         prot_name=prot_df.loc[i, 'prot'],
    #         prot_chain=prot_df.loc[i, 'chain'],
    #         pdb_fp=to('data/pdb/pdbtm/'),
    #         mode='chain',
    #     )
    #     print(threedi_dict)


    prot_df = pd.DataFrame({
        'prot': ['3pux', '3rko', '3udc', '3vr8', '4kjs', '4pi2', ],
        'chain': ['G', 'A', 'A', 'D', 'A', 'C', ],
    })
    for i in prot_df.index:
        print('No.{}: protein: {} chain: {}'.format(i + 1, prot_df.loc[i, 'prot'], prot_df.loc[i, 'chain']))
        # dssp_rsa_run(
        #     prot_name=prot_df.loc[i, 'prot'],
        #     prot_chain=prot_df.loc[i, 'chain'],
        #     pdb_fp='data/pdb/pdbtm/',
        #     sv_fp='data/rsa/',
        # )
        df_rsa = dssp_rsa_access(
            prot_name=prot_df.loc[i, 'prot'],
            prot_chain=prot_df.loc[i, 'chain'],
            rsa_fp=to('data/rsa/')
        )
        print(df_rsa)