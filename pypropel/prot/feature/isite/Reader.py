__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np
import pandas as pd
from pypropel.util.Reader import Reader as pfreader


class Reader:

    def __init__(self, ):
        self.pfreader = pfreader()

    def graphppis(
            self,
            graphppis_fpn : str,
    ):
        df = self.pfreader.generic(
            graphppis_fpn,
            df_sep='\t',
            skiprows=0,
            header=1,
        )
        df = df.reset_index()
        df = df.rename(
            columns={
                'index': 'index',
                'AA': 'aa',
                'Prob': 'pred_prob',
                'Pred': 'pred_label',
            }
        )
        return df

    def mbpred(self, mbp_path, file_name, file_chain, dp=None, by=None, dist_df=None, pair_list=None, sort_=0, is_sort=False, id=0, L=50, len_seq=50):
        self.__sort_ = sort_
        results = self.pfreader.generic(
            mbp_path + file_name + file_chain + '.mbpred',
            df_sep=',',
            header=0,
            is_utf8=True
        )
        results.columns = [
            'index',
            'aa',
            'interact_id',
            'score'
        ]
        # print(results)
        results['aa'] = results['aa'].astype(str)
        recombine = results[[
            'interact_id',
            'score'
        ]]
        # print(recombine)
        # print(recombine.dtypes)
        if self.__sort_ == 1:
            pair_df = pd.DataFrame(pair_list)
            # print(pair_df)
            recombine = self.sort_3(
                recombine,
                is_uniform=True,
                uniform_df=pair_df
            )
            # print(recombine)
            dist_mbpred, dist_true = self.sort_1(recombine, dp, dist_df)
            return dist_mbpred, dist_true
        elif self.__sort_ == 2:
            return self.sort_2(recombine, by=by)
        elif self.__sort_ == 3:
            pass
        elif self.__sort_ == 4:
            pass
        elif self.__sort_ == 5:
            pass
        else:
            return recombine

    def delphi(self, delphi_path, file_name, file_chain, dp=None, by=None, dist_df=None, pair_list=None, sort_=0, is_sort=False, id=0, L=50, len_seq=50):
        self.__sort_ = sort_
        delphi_fpn = delphi_path + file_name + file_chain + '.txt'
        with open(delphi_fpn) as file:
            cues = []
            for line in file:
                if line.split()[0] == '#':
                    continue
                else:
                    cues.append(line.split())
            # print(cues)
            results = pd.DataFrame(cues)
            # print(results)
        results.columns = [
            'interact_id',
            'aa',
            'score'
        ]
        # print(results)
        results['aa'] = results['aa'].astype(str)
        results['interact_id'] = results['interact_id'].astype(np.int)
        results['score'] = results['score'].astype(np.float)
        recombine = results[[
            'interact_id',
            'score'
        ]]
        # print(recombine)
        # print(recombine.dtypes)
        if self.__sort_ == 1:
            pair_df = pd.DataFrame(pair_list)
            # print(pair_df)
            recombine = self.sort_3(
                recombine,
                is_uniform=True,
                uniform_df=pair_df
            )
            # print(recombine)
            dist_delphi, dist_true = self.sort_1(recombine, dp, dist_df)
            return dist_delphi, dist_true
        elif self.__sort_ == 2:
            return self.sort_2(recombine, by=by)
        elif self.__sort_ == 3:
            pass
        elif self.__sort_ == 4:
            pass
        elif self.__sort_ == 5:
            pass
        else:
            return recombine

    def tma300(self,  tma300_path, file_name, file_chain, dp=None, by=None, dist_df=None, pair_list=None, sort_=0, is_sort=False, id=0, L=50, len_seq=50):
        self.__sort_ = sort_
        results = self.pfreader.generic(
            tma300_path + file_name + file_chain + '.tma300',
            df_sep='\t',
            header=None,
            is_utf8=True
        )
        # print(results)
        results.columns = [
            'interact_id',
            'aa',
            'score'
        ]
        results['aa'] = results['aa'].astype(str)
        recombine = results[[
            'interact_id',
            'score'
        ]]
        # print(recombine)
        # print(recombine.dtypes)
        if self.__sort_ == 1:
            pair_df = pd.DataFrame(pair_list)
            # print(pair_df)
            recombine = self.sort_3(
                recombine,
                is_uniform=True,
                uniform_df=pair_df
            )
            # print(recombine)
            dist_tma300, dist_true = self.sort_1(recombine, dp, dist_df)
            return dist_tma300, dist_true
        elif self.__sort_ == 2:
            return self.sort_2(recombine, by=by)
        elif self.__sort_ == 3:
            pass
        elif self.__sort_ == 4:
            pass
        elif self.__sort_ == 5:
            pass
        else:
            return recombine


if __name__ == "__main__":
    from pypropel.path import to

    p = Reader()

    print(p.graphppis(
        graphppis_fpn=to('')
    ))

