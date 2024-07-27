__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import pandas as pd


class Reader:

    def generic(
            self,
            df_fpn,
            df_sep='\t',
            skiprows=None,
            header=None,
            encoding='utf-8',
            comment=None
    ):
        return pd.read_csv(
            df_fpn,
            sep=df_sep,
            header=header,
            encoding=encoding,
            skiprows=skiprows,
            comment=comment,
        )


    def excel(self, df_fpn, sheet_name='Sheet1', header=None, is_utf8=False):
        if is_utf8:
            return pd.read_excel(
                df_fpn,
                sheet_name=sheet_name,
                header=header,
                encoding='utf-8',
                engine='openpyxl',
            )
        else:
            return pd.read_excel(
                df_fpn,
                sheet_name=sheet_name,
                header=header,
                engine='openpyxl',
            )