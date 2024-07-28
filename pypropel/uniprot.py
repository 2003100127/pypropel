__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

import pandas as pd
from pypropel.prot.uniprot.Text import Text as uniprottext


def from_text(
        text_fpn : str,
        sv_json_fpn : str,
        sv_df_fpn : str,
) -> pd.DataFrame:
    return uniprottext(
        text_fpn=text_fpn,
        sv_json_fpn=sv_json_fpn,
        sv_df_fpn=sv_df_fpn,
    ).parse()


if __name__ == "__main__":
    from pypropel.path import to

    print(from_text(
        text_fpn=to('data/uniprot/text/uniprotkb_Human_AND_reviewed_true_AND_m_2023_11_29.txt'),
        sv_json_fpn=to('data/uniprot/text/human.json'),
        sv_df_fpn=to('data/uniprot/text/human.txt'),
    ))