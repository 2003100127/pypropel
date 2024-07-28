__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import linecache
import pandas as pd
from pypropel.util.Reader import Reader as pfreader


class Reader:

    def __init__(self, ):
        self.pfreader = pfreader()

    def spider3(
            self,
            spider3_path,
            prot_name,
            file_chain,
    ):
        spider3 = self.pfreader.generic(
            df_fpn=spider3_path + prot_name + file_chain + '.spd33',
            df_sep='\s+',
            header=0,
        )
        # print(spider3.dtypes)
        spider3[[
            'ASA',
            'Phi',
            'Psi',
            'Theta(i-1=>i+1)',
            'Tau(i-2=>i+2)',
            'HSE_alpha_up',
            'HSE_alpha_down',
            'HSE_beta_up',
            'HSE_beta_down',
            'CN',
        ]] = spider3[[
            'ASA',
            'Phi',
            'Psi',
            'Theta(i-1=>i+1)',
            'Tau(i-2=>i+2)',
            'HSE_alpha_up',
            'HSE_alpha_down',
            'HSE_beta_up',
            'HSE_beta_down',
            'CN',
        ]].apply(lambda x: x / 100)
        return spider3

    def spider3_to_ss(self, spider3_path, prot_name, file_chain, sv_fp, state=3):
        ss_aa_spider3 = self.spider3(
            spider3_path=spider3_path,
            prot_name=prot_name,
            file_chain=file_chain
        )['SS']
        if state == 3:
            with open(sv_fp + str(prot_name) + str(file_chain) + '.ss', 'w') as f:
                f.write('>' + str(prot_name) + str(file_chain) + '\n')
                f.write(''.join(str(x) for x in ss_aa_spider3.values.tolist()))
                f.close()
            return ss_aa_spider3

    def psipred(self, prot_name, file_chain, psipred_ss_path=None, psipred_ss2_path=None, psipred_horiz_path=None, ):
        if psipred_ss_path is not None:
            psipred_ss = self.pfreader.generic(
                df_fpn=psipred_ss_path + prot_name + file_chain + '.ss',
                df_sep='\s+',
                header=None,
            )
            return psipred_ss
        if psipred_ss2_path is not None:
            f = open(psipred_ss2_path + prot_name + file_chain + '.ss2', 'r')
            psipred_ss2 = pd.DataFrame()
            c = 0
            for line in f.readlines():
                line = line.split()
                if len(line) and line[0] is not '#':
                    psipred_ss2[c] = line
                    c += 1
            psipred_ss2_ = psipred_ss2.T
            psipred_ss2_[0] = psipred_ss2_[0].astype(int)
            psipred_ss2_[[3, 4, 5]] = psipred_ss2_[[3, 4, 5]].astype(float)
            return psipred_ss2_
        if psipred_horiz_path is not None:
            psipred_horiz = self.pfreader.generic(
                df_fpn=psipred_horiz_path + prot_name + file_chain + '.horiz',
                df_sep='\s+',
                header=0,
            )
            return psipred_horiz

    def sspro(self, sspro_path, prot_name, file_chain):
        line = linecache.getline(sspro_path + prot_name + file_chain + '.ss', 2)
        sspro_ = pd.DataFrame(list(line)[:-1])
        return sspro_

    def sspro8(self, sspro8_path, prot_name, file_chain):
        line = linecache.getline(sspro8_path + prot_name + file_chain + '.ss8', 2)
        sspro8_ = pd.DataFrame(list(line)[: -1])
        return sspro8_


if __name__ == "__main__":
    from pypropel.path import to

    p = Reader()

    # print(p.spider3(
    #     spider3_path=to('data/ss/spider3/'),
    #     prot_name='E',
    #     file_chain=''
    # ))

    # print(p.spider3_to_ss(
    #     spider3_path=to('data/ss/spider3/'),
    #     prot_name='E',
    #     file_chain='',
    #     sv_fp=to('data/ss/spider3/'),
    # ))

    # print(p.psipred(
    #     # psipred_ss_path=to('data/ss/psipred/'),
    #     # psipred_ss2_path=to('data/ss/psipred/'),
    #     psipred_horiz_path=to('data/ss/psipred/'),
    #     prot_name='1aig', # 1aig L 1bcc D
    #     file_chain='L',
    # ))

    # print(p.sspro(
    #     sspro_path=to('data/ss/sspro/'),
    #     prot_name='1aig',
    #     file_chain='L'
    # ))

    # print(p.sspro8(
    #     sspro8_path=to('data/ss/sspro8/'),
    #     prot_name='1aig',
    #     file_chain='L'
    # ))
