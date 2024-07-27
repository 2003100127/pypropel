__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"


class Standardize:

    def __init__(self, ):
        pass

    def minmax(self, value, min, max):
        return (value - min) / max - min

    def minmax1(self, arr):
        min_ = min(arr)
        max_ = max(arr)
        return [(i-min_)/(max_-min_) for i in arr]

    def minmax2(self, idict):
        max_id = max(idict, key=idict.get)
        min_id = min(idict, key=idict.get)
        # print(max_id)
        # print(min_id)
        max_ = idict[max_id]
        min_ = idict[min_id]
        for k, v in idict.items():
            idict[k] = (v - min_) / (max_ - min_)
        return idict


if __name__ == "__main__":
    p = Standardize()
    idd = {
        'A': 31, 'C': 55, 'D': 54, 'E': 83,
        'F': 132, 'G': 3, 'H': 96, 'I': 111,
        'K': 119, 'L': 111, 'M': 105, 'N': 56,
        'P': 32.5, 'Q': 85, 'R': 124, 'S': 32,
        'T': 61, 'V': 84, 'W': 170, 'Y': 136,
    }
    print(p.minmax2(idd))