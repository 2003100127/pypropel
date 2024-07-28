__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from Bio import AlignIO
from pypropel.util.Reader import Reader as pfreader
from pypropel.util.Console import Console


class Convert:

    def __init__(
            self,
            input_fpn,
            output_fpn: str = None,
            output_fp : str = None,
            in_format='fasta',
            out_format: str = 'clustal',
            verbose: bool = True,
    ):
        self.input_fpn = input_fpn
        self.output_fpn = output_fpn
        self.in_format = in_format
        self.out_format = out_format
        self.output_fp = output_fp

        self.pfreader = pfreader()
        self.console = Console()
        self.console.verbose = verbose

    def read_msa(self, ):
        return AlignIO.read(handle=self.input_fpn, format=self.in_format)

    def reformat(
            self,
    ):
        """

        Attributes
        ----------
        in_format or out_format that it supports
            clustal
            emboss
            fasta
            fasta-m10
            ig
            msf
            nexus
            phylip
            phylip-sequential
            phylip-relaxed
            stockholm
            mauve

        Returns
        -------

        """
        self.console.print("=========> convert from {} to {}".format(self.in_format, self.out_format))
        AlignIO.convert(
            self.input_fpn,
            self.in_format,
            out_file=self.output_fpn,
            out_format=self.out_format,
        )
        return 'Finished'

    def tofasta_sgl(
            self,
    ) -> str:
        f = open(self.input_fpn)
        msa = [[str(line) for line in f_parser.split()] for f_parser in f]
        ids = msa[::2]
        homologs = msa[1::2]
        # import os
        # output_fp = os.path.dirname(os.path.realpath(self.output_fpn))
        for id, homolog in zip(ids, homologs):
            self.console.print("=========>extract {} to save".format(id[0][1:]))
            with open(self.output_fp + '/' + id[0][1:] + '.fasta', 'w') as f_renew:
                # print(self.output_fp + id[0][1:] + '.fasta')
                f_renew.write(id[0] + '\n')
                f_renew.write(homolog[0] + '\n')
        return 'Finished'

    def write(self, input, output_fpn, out_format='fasta'):
        with open(output_fpn, "w") as handle:
            AlignIO.write(input, handle, out_format)
        return


if __name__ == "__main__":
    from pypropel.path import to

    p = Convert(
        input_fpn=to('data/msa/a2m/ET.a2m'),
        in_format='fasta',
        output_fpn=to('data/msa/a2m/ET_converted.clustal'),
        out_format='clustal',
    )

    # print(p.reformat())
    print(p.tofasta_sgl())