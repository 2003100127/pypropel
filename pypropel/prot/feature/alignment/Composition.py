__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np
from pypropel.prot.feature.alignment.frequency.Single import Single as fs
# from protein.tool.alignment.Reader import reader as pareader
from pypropel.util.Console import Console


class Composition:

    def __init__(
            self,
            msa,
            verbose: bool = True,
    ):
        """
        Amino acid composition for a sequence with its MSA is calculated.

        Parameters
        ----------
        msa
        """
        self.msa = msa
        # self.pareader = pareader()
        self.fs = fs(self.msa)

        self.console = Console()
        self.console.verbose = verbose

    def aac(self, ):
        """
        It gets amino acid composition SingleFrequency

        Returns
        -------
            Dict

        """
        return self.fs.alignment()

    def evolutionary_profile(self, ):
        """
        Notes
        -----
        Evolutionary profile is a term of position-specific
        scoring matrix (PSSM), which is calculated with the
        most general way. Evolutionary profile is calculated by
        sigmoid(singleFreq2Col() / singleFreq2MSA()).

        Examples
        --------
        1atzA: QPLDVILLLDGSSSFPASYFDEMKSFAKAFISKANIGPRLTQVSVLQYGSITTIDVPWNVVPEKAHLLSLVDVMQ
            :return 75x21 matrix

        References
        ----------
        Hönigschmid, P., & Frishman, D. (2016). Accurate prediction of helix
            interactions and residue contacts in membrane protein. Journal of
            structural biology, 194(1), 112-123.

        Returns
        -------
            numpy.ndarray
            2d array - row: length of sequence; col: 21

        """
        freq_column, _ = self.fs.columns()
        # print(freq_column)
        aac = self.fs.alignment()
        fraction = {}
        for aa in self.fs.aa_alphabet:
            if aac[aa] == 0:
                fraction[aa] = freq_column[aa] * 0
            else:
                fraction[aa] = freq_column[aa] / aac[aa]
        return fraction

    def evolutionary_profile_norm(self, ):
        fraction = self.evolutionary_profile()
        for aa in self.fs.aa_alphabet:
            for j in range(len(fraction['A'])):
                if fraction[aa][j] == 0:
                    fraction[aa][j] = 1
                else:
                    part = -np.log(fraction[aa][j])
                    fraction[aa][j] = 1 / (1 + np.exp(part))
        return fraction




if __name__ == "__main__":
    from pypropel.prot.feature.alignment.MSA import MSA as msaparser
    from pypropel.path import to

    msa = msaparser(
        msa_fpn=to('data/msa/aln/1aijL.aln'),
    ).read()

    p = Composition(
        msa=msa,
    )

    # print(p.aac())

    # print(p.evolutionary_profile())

    # print(len(p.evolutionary_profile_norm()['A']))

