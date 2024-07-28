__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import pandas as pd
from pypropel.util.Reader import Reader as pfreader
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Console import Console


class ListDiffer:

    def __init__(
            self,
            verbose: bool = True,
    ):
        self.pfreader = pfreader()
        self.pfwriter = pfwriter()
        self.console = Console()
        self.console.verbose = verbose

    def unipartite(
            self,
            pds_lg,
            pds_sm,
            sv_diff_fpn=None,
            sv_rept_fpn=None,
    ):
        """
        Examples
        --------
        # slow mode
        # for i in range(df_sm.shape[0]):
        #     query = df_lg.loc[df_lg[col] == df_sm[col][i]]
        #     if len(query.index):
        #         row = query.index[0]
        #         repeated.append(df_sm[col][i])
        #     else:
        #         continue
        #     df_lg.drop([row], inplace=True)

        Parameters
        ----------
        pds_lg
            the 1st longer pandas Series
        pds_sm
            the 2nd shorter pandas Series
        sv_diff_fpn
            file contains different content between two lists
        sv_rept_fpn
            file contains repeated content between two lists

        Returns
        -------

        """
        lg_set = set(pds_lg.values.tolist())
        sm_set = set(pds_sm.values.tolist())
        pds_differ = pd.Series(list(lg_set.difference(sm_set)))
        pds_repeat = pd.Series(list(lg_set.intersection(sm_set)))
        if sv_diff_fpn:
            self.pfwriter.generic(pds_differ, sv_diff_fpn, df_sep='\t')
        if sv_rept_fpn:
            self.pfwriter.generic(pds_repeat, sv_rept_fpn, df_sep='\t')
        return pds_differ, pds_repeat

    def bipartite(
            self,
            pds_lg_1,
            pds_lg_2,
            pds_sm_1,
            pds_sm_2,
            sv_diff_fpn=None,
            sv_rept_fpn=None,
    ):
        """
        df_lg_agent = df_lg
        # slow mode 1
            for j, _ in df_sm.iterrows():
                for i, _ in df_lg.iterrows():
                    if df_sm.loc[j][0] + df_sm.loc[j][1] == df_lg.loc[i][0] + df_lg.loc[i][1]:
                        df_lg_agent.drop([i], inplace=True)
            df_lg_agent = df_lg_agent.reset_index(drop=True)
        # slow mode 2
            same = []
            for j, _ in df_sm.iterrows():
                for i, _ in df_lg.iterrows():
                    if df_sm.loc[j][0]+df_sm.loc[j][1] == df_lg.loc[i][0]+df_lg.loc[i][1]:
                        same.append([df_sm.loc[j][0] + ' ' + df_sm.loc[j][1]])
            same = pd.DataFrame(same)
        # slow mode 3 (fast a bit)
            repeated = []
            for i in range(df_sm.shape[0]):
                print(i)
                query = df_lg.loc[
                    (df_lg[0] == df_sm[0][i]) & (df_lg[1] == df_sm[1][i])
                ]
                if len(query.index):
                    row = query.index[0]
                    repeated.append([df_sm[0][i], df_sm[1][i]])
                else:
                    continue
                df_lg_agent.drop([row], inplace=True)
            df_lg_agent = df_lg_agent.reset_index(drop=True)
            if sv_diff_fpn:
                self.pfwriter.generic(df_lg_agent, sv_diff_fpn, df_sep='\t')
            if sv_rept_fpn:
                self.pfwriter.generic(repeated, sv_rept_fpn, df_sep='\t')

        Parameters
        ----------
        pds_lg_1
            the 1st longer pandas Series contains the 1st column
        pds_lg_2
            the 2nd shorter pandas Series contains the 2nd column
        pds_sm_1
            the 1st longer pandas Series contains the 1st column
        pds_sm_2
            the 2nd shorter pandas Series contains the 2nd column
        sv_diff_fpn
            file contains different content between two lists
        sv_rept_fpn
            file contains repeated content between two lists

        Returns
        -------

        """
        df_lg = pds_lg_1.to_frame().join(pds_lg_2)
        df_sm = pds_sm_1.to_frame().join(pds_sm_2)
        pds_differ, pds_repeat = self.unipartite(
            pds_lg=df_lg.apply(lambda x: str(x[0]) + "_" + str(x[1]), axis=1),
            pds_sm=df_sm.apply(lambda x: str(x[0]) + "_" + str(x[1]), axis=1),
            sv_diff_fpn=sv_diff_fpn,
            sv_rept_fpn=sv_rept_fpn,
        )
        df_differ = pd.DataFrame()
        df_repeat = pd.DataFrame()
        if not pds_differ.empty:
            df_differ[0], df_differ[1] = zip(*pds_differ.apply(lambda x: (x.split("_")[0], x.split("_")[1])))
        if not pds_repeat.empty:
            df_repeat[0], df_repeat[1] = zip(*pds_repeat.apply(lambda x: (x.split("_")[0], x.split("_")[1])))
        return df_differ, df_repeat