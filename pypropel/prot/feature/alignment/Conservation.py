__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import re
import numpy as np
import pandas as pd
from pypropel.util.Reader import Reader as pfreader


class Conservation:
    
    def __init__(self, ):
        self.pfreader = pfreader()

    def get(self, ent_dict, thres=20):
        conser = {}
        for k, v in ent_dict.items():
            conser[k] = 1 - (v / np.log(thres))
        return conser

    def getFromOutside(self, entropy, thres=20):
        # entropy, _ = self.it.entropyMirror(gap_thres=0.5)
        conser = {}
        for k, v in entropy.items():
            conser[k] = 1 - (v / np.log(thres))
        return conser

    def jsd(self, jsd_fpn):
        """
        Jensen–Shannon (JS) divergence

        Examples
        --------
        An JSD file looks like
        # ./CLEC2B_LOC113845378.clustal -- js_divergence - window_size: 3 - window lambda: 0.50 - background: blosum62 - seq. weighting: True - gap penalty: 1 - normalized: False
        # align_column_number	score	column
        0	-1000.000000	M-M-M-T----TET--TTTMM-M-----------------------------------------------------P-----------------------------------------------
        1	-1000.000000	E-S-D-Q----QTE--QQQEE-E-----------------------------------------------------V----------------------G------------------------
        2	-1000.000000	K-S-S-N----DGY--DNNPP-P-----------------------------------------------------P----------------------P------------------------
        3	-1000.000000
        ...
        205	-1000.000000	PN---S---P--VPF------R---------------L---LL---LLL------L----DD---V----A--E----------------------------------EME---L-EE-I----
        206	-1000.000000	EP---P---R---VP----------------------S---SS---SSS------S----HH---T----E--K----------------------------------SCS---C-SS-M----

        References
        ----------
        John A. Capra, Mona Singh, Predicting functionally important residues from
        sequence conservation, Bioinformatics, Volume 23, Issue 15, August 2007,
        Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270


        Parameters
        ----------
        jsd_fpn
            JSD file

        Returns
        -------

        """
        df = self.pfreader.generic(jsd_fpn, df_sep='\s+', comment='#')
        df = df.rename(
            columns={
                0: 'alignment_col',
                1: 'score',
                2: 'seq',
            }
        )
        df['score'] = df['score'].apply(lambda x: 0 if x == -1000.0 else x)
        return df

    def consurf_v1(
            self,
            consurf_fpn : str,
    ):
        """
        v1 means the output adopted by ConSurf before 2024.
        """
        f = open(consurf_fpn, 'r')
        c = []
        for line in f.readlines():
            line = line.strip().split('\t')
            if len(line) != 0 and re.match(r'[0-9]', line[0]):
                c.append([i.strip() for i in line])
        f.close()
        l = pd.DataFrame(c)
        print(l)
        l[0] = l[0].astype(np.int64)
        l[2] = l[2].astype(np.float64)
        l[2] = l[2]
        l = l[[0, 1, 2, 4, 11, 12]]
        l = l.rename(columns={
            0: 'position',
            1: 'amino acid',
            2: 'score',
            4: 'color',
            11: 'exposed/buried',
            12: 'structral/functional',
        })
        return l


if __name__ == "__main__":
    from pypropel.path import to

    p = Conservation()

    print(p.jsd(jsd_fpn=to('data/conservation/jsd/SR24_CtoU/CLEC2B_LOC113845378.jsd')))

    # print(p.consurf_v1(
    #     to('data/conservation/consurf/E.consurf')
    # ))