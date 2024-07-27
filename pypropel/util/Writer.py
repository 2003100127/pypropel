__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import pandas as pd


class Writer:

    def generic(self, df, sv_fpn, df_sep='\t', header=None, index=False, id_from=0):
        df_ = pd.DataFrame(df)
        df_.index = df_.index + id_from
        return df_.to_csv(
            sv_fpn,
            sep=df_sep,
            header=header,
            index=index
        )

    def excel(self, df, sv_fpn=None, sheet_name='Sheet1', header=None, index=False, id_from=0):
        df_ = pd.DataFrame(df)
        df_.index = df_.index + id_from
        return df_.to_excel(
            sv_fpn,
            sheet_name=sheet_name,
            header=header,
            index=index
        )