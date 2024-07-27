__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"


class Symbol:

    def __init__(self, ):
        self.aa = self.single(gap=False)

    def single(self, gap=False, universal=False):
        if universal:
            if gap:
                return ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
            else:
                return ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
        else:
            if gap:
                return ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '-']
            else:
                return ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

    def pair(self, ):
        paris_ = []
        for _, i in enumerate(self.aa):
            for _, j in enumerate(self.aa):
                    paris_.append(i + j)
        return paris_

    def triple(self, ):
        kmers_ = []
        for _, i in enumerate(self.aa):
            for _, j in enumerate(self.aa):
                for _, k in enumerate(self.aa):
                    kmers_.append(i + j + k)
        return kmers_

    def quadruple(self, ):
        kmers_ = []
        for _, i in enumerate(self.aa):
            for _, j in enumerate(self.aa):
                for _, k in enumerate(self.aa):
                    for _, l in enumerate(self.aa):
                        kmers_.append(i + j + k + l)
        return kmers_

    def todict(self, gap=False):
        aa_dict = {}
        for k, v in enumerate(self.single(gap=gap)):
            aa_dict[v] = k
        return aa_dict


if __name__ == "__main__":
    p = Symbol()
    # print(p.single())
    # print(p.pair())
    # print(p.triple())