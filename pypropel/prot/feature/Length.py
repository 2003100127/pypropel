__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from typing import List

import numpy as np


class Length:

    def __init__(self, ):
        pass

    def segment(
            self,
            seg_lower : List,
            seg_upper : List,
    ):
        snippets = []
        for i in range(len(seg_lower)):
            snippets.append(seg_upper[i] - seg_lower[i] + 1)
        len_seg = sum(snippets)
        return len_seg

    def msa(
            self,
            msa,
    ):
        return len(msa)

    def sequence(
            self,
            sequence,
    ):
        return len(sequence)

    def log(self, sequence):
        return np.log(self.sequence(sequence))