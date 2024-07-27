__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np
from scipy import stats
from sklearn.preprocessing import normalize


class Activate:

    def __init__(self, ):
        pass
    
    def sigmoid(self, value):
        return 1 / (1 + np.exp(-value))

    def sigmoid1(self, arr):
        pass

    def sigmoid2(self, idict):
        for k, v in idict.items():
            idict[k] = self.sigmoid(v)
        return idict

    def zscore(self, idict):
        idict_ = {}
        keys = idict.keys()
        zscores = stats.zscore(list(idict.values()))
        for i, key in enumerate(list(keys)):
            idict_[key] = zscores[i]
        return idict_

    def sklearn(self, idict):
        idict_ = {}
        keys = idict.keys()
        x = np.array(list(idict.values()))
        sklearn = normalize(x[:, np.newaxis], axis=0).ravel()
        for i, key in enumerate(list(keys)):
            idict_[key] = sklearn[i]
        return idict_


if __name__ == "__main__":
    p = Activate()
    idd = {
        'A': 31, 'C': 55, 'D': 54, 'E': 83,
        'F': 132, 'G': 3, 'H': 96, 'I': 111,
        'K': 119, 'L': 111, 'M': 105, 'N': 56,
        'P': 32.5, 'Q': -85, 'R': 124, 'S': 32,
        'T': 61, 'V': 84, 'W': 170, 'Y': 136
    }
    print(p.sigmoid(62.92868827638704))
    print(p.sigmoid2(idd))