__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from Bio import SeqIO


class Fasta:

    def get(
            self,
            fasta_fpn: str,
    ) -> str:
        """
        Get the sequence from a FASTA file.

        Parameters
        ----------
        fasta_fpn : str
            The filepath of the FASTA file.

        Returns
        -------
        str
            A sequence from the FASTA file.

        Raises
        ------
        EmptyFastaError
            If the FASTA file does not contain any sequence.
        """
        sequence = "".join([str(seq.seq) for seq in SeqIO.parse(fasta_fpn, "fasta")])
        return sequence

    def save(
            self,
            list_2d,
            sv_fp,
    ):
        for i, e in enumerate(list_2d):
            prot_name = str(e[0])
            seq = str(e[1])
            print('No.{} saving {} in FASTA format.'.format(i+1, prot_name))
            f = open(sv_fp + prot_name + '.fasta', 'w')
            f.write('>' + prot_name + '\n')
            f.write(seq + '\n')
            f.close()
        return 0