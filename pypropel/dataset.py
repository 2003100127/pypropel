__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import Tuple

import pandas as pd

from pypropel.prot.file.Pack import Pack


def download_pack(
        prot_df : pd.DataFrame,
        pdb_cplx_fp : str,
        pdb_fp : str,
        xml_fp : str,
        fasta_fp : str,
):
    return Pack(
        prot_df=prot_df,
    ).execute(
        pdb_cplx_fp=pdb_cplx_fp,
        pdb_fp=pdb_fp,
        xml_fp=xml_fp,
        fasta_fp=fasta_fp,
    )


if __name__ == "__main__":
    from pypropel.path import to

    import pandas as pd

    prot_df = pd.DataFrame({
        'prot': ['1aig', '1aij', '1xqf'],
        'chain': ['L', 'L', 'A'],
    })

    print(download_pack(
        prot_df=prot_df,
        pdb_cplx_fp=to('data/pdb/complex/pdbtm/'),
        pdb_fp=to('data/tmp/'),
        xml_fp=to('data/tmp/'),
        fasta_fp=to('data/tmp/'),
    ))