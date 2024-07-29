__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"


class MSA:

    def __init__(
            self,
            msa_fpn,
    ):
        """
        It provides elementary operations of parsing a MSA

        Parameters
        ----------
        msa_fpn
            msa's file path and file name

        """
        self.msa_fpn = msa_fpn

    def read(self, ) -> list:
        """
        parse a MSA from an aln file.

        Examples
        --------
           a fasta file
           >sequence|1atz|chain:A
           QPLDVILLLDGSSSFPASYFDEMKSFAKAFISKANIGPRLTQVSVLQYGSITTIDVPWNVVPEKAHLLSLVDVMQ

           return
            QPLDVILLLDGSSSFPASYFDEMKSFAKAFISKANIGPRLTQVSVLQYGSITTIDVPWNVVPEKAHLLSLVDVMQ

        Returns
        -------
        1d array : np ndarray
            row for number of multiple sequences; column: length of protein sequence.

        """
        read_msa = open(self.msa_fpn, 'r')
        results = list()
        for line in read_msa.readlines():
            line = line.strip()
            if not len(line) or line.startswith('#'):
                continue
            results.append(line)
        return results