__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2020"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"


class Symbol:

    def __init__(self, ):
        self.nt = self.single(gap=False)

    def single(self, gap=False, universal=False):
        if universal:
            if gap:
                return []
            else:
                return []
        else:
            if gap:
                return ['A', 'C', 'G', 'U', '-']
                # return ['A', 'C', 'G', 'T', '-']
            else:
                return ['A', 'C', 'G', 'U']
                # return ['A', 'C', 'G', 'T']

    def pair(self, ):
        paris_ = []
        for _, i in enumerate(self.nt):
            for _, j in enumerate(self.nt):
                    paris_.append(i + j)
        return paris_

    def triple(self, ):
        kmers_ = []
        for _, i in enumerate(self.nt):
            for _, j in enumerate(self.nt):
                for _, k in enumerate(self.nt):
                    kmers_.append(i + j + k)
        return kmers_

    def quadruple(self, ):
        kmers_ = []
        for _, i in enumerate(self.nt):
            for _, j in enumerate(self.nt):
                for _, k in enumerate(self.nt):
                    for _, l in enumerate(self.nt):
                        kmers_.append(i + j + k + l)
        return kmers_

    def todict(self, gap=False):
        nt_dict = {}
        for k, v in enumerate(self.single(gap=gap)):
            nt_dict[v] = k
        return nt_dict


if __name__ == "__main__":
    p = Single()
    print(p.single())
    print(p.todict())