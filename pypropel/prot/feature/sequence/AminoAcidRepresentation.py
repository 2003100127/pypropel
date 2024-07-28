__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from typing import  List

import time
from pypropel.prot.sequence.Symbol import Symbol
from pypropel.util.Console import Console


class AminoAcidRepresentation:
    
    def __init__(
            self,
            verbose : bool = True,
    ):
        self.aa_dict = Symbol().todict(gap=False)
        self.console = Console()
        self.console.verbose = verbose
        
    def onehot(
            self,
            arr_2d : List[List],
            arr_aa_names : List[List],
    ) -> List[List]:
        start_time = time.time()
        for i, aa_win_names in enumerate(arr_aa_names):
            for j in aa_win_names:
                if j is None:
                    for k in range(20):
                        arr_2d[i].append(0)
                else:
                    bool_ = [0] * 20
                    bool_[self.aa_dict[j]] = 1
                    for k in range(20):
                        arr_2d[i].append(bool_[k])
        end_time = time.time()
        self.console.print("=========>AA representation: {time}s.".format(time=end_time - start_time))
        return arr_2d