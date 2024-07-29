__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import re
import linecache
import pandas as pd
from pypropel.util.Reader import Reader as pfreader


class Reader:

    def __init__(self, ):
        self.pfreader = pfreader()

    def solvpred(
            self,
            solvpred_fp,
            prot_name,
            file_chain,
    ):
        df = self.pfreader.generic(
            df_fpn=solvpred_fp + prot_name + file_chain + '.solv',
            df_sep='\s+',
            header=None
        )
        return df

    def accpro(
            self,
            accpro_fp,
            prot_name,
            file_chain,
    ):
        line = linecache.getline(accpro_fp + prot_name + file_chain + '.acc', 2)
        df = pd.DataFrame(list(line)[:-1])
        return df

    def accpro20(
            self,
            accpro20_fp,
            prot_name,
            file_chain,
    ):
        line = linecache.getline(accpro20_fp + prot_name + file_chain + '.acc20', 2)
        df = pd.DataFrame([float(e) for e in re.split(r' ', line)])
        # print(accpro20)
        return df.apply(lambda x: x/100)


if __name__ == "__main__":
    from pypropel.path import to

    p = Reader()

    print(p.solvpred(
        solvpred_fp=to('data/accessibility/solvpred/'),
        prot_name='1aig',
        file_chain='L',
    ))

    print(p.accpro(
        accpro_fp=to('data/accessibility/accpro/'),
        prot_name='1aig',
        file_chain='L',
    ))

    print(p.accpro20(
        accpro20_fp=to('data/accessibility/accpro20/'),
        prot_name='1aig',
        file_chain='L',
    ))