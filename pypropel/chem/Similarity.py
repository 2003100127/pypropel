__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"

import pandas as pd
from rdkit import DataStructs
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Console import Console


class Similarity:

    def __init__(self, verbose=False):
        self.console = Console()
        self.verbose = verbose
        self.console.verbose = self.verbose

    def get(self, fps, type='Tanimoto', is_sv=False, sv_fpn='./'):
        """

        :param fps: a list of fingerprints - self.fingerprints(type='Morgan')
        :param type: Tanimoto, Dice, Cosine, Sokal, Russel,
                     Kulczynski, McConnaughey, and Tversky

        :return:
        """
        self.console.print('===>calculating similarity using {}'.format(type))
        num_fps = len(fps)
        if type == 'Tanimoto':
            method = DataStructs.TanimotoSimilarity
        elif type == 'Dice':
            method = DataStructs.DiceSimilarity
        elif type == 'Sokal':
            method = DataStructs.SokalSimilarity
        elif type == 'Russel':
            method = DataStructs.RusselSimilarity
        elif type == 'Kulczynski':
            method = DataStructs.KulczynskiSimilarity
        elif type == 'McConnaughey':
            method = DataStructs.McConnaugheySimilarity
        else:
            method = DataStructs.TverskySimilarity
        sims = []
        for i in range(num_fps):
            if i!= 0 and i % 1000 == 0:
                self.console.print('======>{} molecules processed'.format(i))
            for j in range(i+1, num_fps):
                sims.append([
                    fps[i][0],
                    fps[j][0],
                    DataStructs.FingerprintSimilarity(
                        fps[i][1], fps[j][1], metric=method
                    ),
                ])
        sims_ = pd.DataFrame(sims)
        if is_sv:
            pfwriter().generic(sims_, sv_fpn)
        return sims_


if __name__ == "__main__":
    from pypropel.path import to
    from pypropel.chem.SDF import SDF
    from pypropel.chem.Fingerprint import Fingerprint

    p1 = Similarity()
    p2 = Fingerprint()
    p3 = SDF(sdf_fpn=to('data/chem/approved.sdf'))
    # p3.kind_ = 'all_2d'
    # p3.kind_ = 'approved_2d'
    p3.kind = 'expt_2d'
    mols = p3.extract()
    print(len(mols))
    fps = p2.get(mols, is_to_vec=False, single=False)
    print(fps)

    print(p1.get(
        fps=fps,
        type='Tanimoto',
        is_sv=True,
        sv_fpn=to('data/chem/expt_2d.si'),
    ))
