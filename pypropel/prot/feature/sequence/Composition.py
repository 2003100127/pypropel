__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from collections import Counter
from pypropel.prot.sequence.Symbol import Symbol
from pypropel.util.Console import Console


class Composition:
    """
    Amino acid composition for a sequence with its MSA is calculated.

    """
    def __init__(
            self,
            sequence,
            verbose: bool = True,
    ):
        self.sequence = sequence
        self.symbol = Symbol()
        self.aa = self.symbol.single()

        self.console = Console()
        self.console.verbose = verbose

    def aac(self, ):
        aac_ = {}
        for _, i in enumerate(self.aa):
            aac_[i] = round(self.sequence.count(i) / len(self.sequence), 6)
        return aac_

    def dac(self, ):
        """
        dipeptide composition (DPC).
        """
        dac_ = []
        for _, i in enumerate(self.symbol.pair()):
            dac_.append([i, round(self.sequence.count(i) / (len(self.sequence) - 1), 6)])
        return dac_

    def tac(self, ):
        tac_ = []
        for _, i in enumerate(self.symbol.triple()):
            tac_.append([i, round(self.sequence.count(i) / (len(self.sequence) - 2), 6)])
        return tac_

    def qac(self, ):
        qac_ = []
        for _, i in enumerate(self.symbol.quadruple()):
            qac_.append([i, round(self.sequence.count(i) / (len(self.sequence) - 3), 6)])
        return qac_

    def cksnap(self, k):
        """
        Composition of k-Spaced amino acid pairs (CKSAAP - 16)

        Parameters
        ----------
        k

        Returns
        -------

        """
        c = []
        for i in range(len(self.sequence) - (k + 1)):
            c.append(self.sequence[i] + self.sequence[i + k + 1])
        freq_dict = dict(Counter(c))
        cksnap_ = {}
        # print(freq_dict)
        keys = [*freq_dict.keys()]
        for _, i in enumerate(self.symbol.pair()):
            if i in keys:
                cksnap_[i] = round(freq_dict[i] / (len(self.sequence) - (k + 1)), 6)
                # cksnap_[i] = freq_dict[i]
            else:
                cksnap_[i] = 0
        return cksnap_

    def aveanf(self, ):
        """average Accumulated Amino Acid Frequency (AAAF)"""
        cdict = {}
        for aa_ref in self.aa:
            cdict[aa_ref] = []
            cnt = 0
            for i, nt in enumerate(self.sequence):
                if nt == aa_ref:
                    cnt += 1
                    cdict[aa_ref].append(cnt / (i + 1))
        # print(cdict)
        cdict_ = {}
        for k, v in cdict.items():
            if v != []:
                cdict_[k] = sum(v) / len(v)
            else:
                cdict_[k] = 0
        return cdict_


if __name__ == "__main__":
    seq = "ADGCGVGEGTGQGPMCNCMCMKWVYADEDAADLESDSFADEDASLESDSFPWSNQRVFCSFADEDAS"

    p = Composition(
        sequence=seq
    )
    # print(p.aac())

    # print(p.dac())

    # print(p.tac())

    # print(p.qac())

    # print(p.cksnap(2))

    # print(p.aveanf())
