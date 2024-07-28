__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from pypropel.prot.sequence.Fasta import Fasta
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Console import Console


class IsEmpty:

    def __init__(
            self,
            prot_df,
            sv_empty_fp,
            fasta_fp,
            verbose: bool = True,
    ):
        self.prot_df = prot_df
        self.fasta_fp = fasta_fp
        self.sv_empty_fp = sv_empty_fp
        self.pfwriter = pfwriter()
        self.console = Console()
        self.console.verbose = verbose

    def fasta(
            self,
    ) -> str:
        empty = []
        for i in range(self.prot_df.shape[0]):
            prot_name = self.prot_df['prot'][i]
            prot_chain = self.prot_df['chain'][i]
            file_chain = chainname().chain(prot_chain)
            self.console.print('============>No{}. protein {} chain {}'.format(i, prot_name, prot_chain))
            try:
                if Fasta().get(fasta_fpn=self.fasta_fp + prot_name + file_chain + '.fasta') == '':
                    self.console.print('===============>The fasta file is empty.')
                    empty.append([prot_name, prot_chain])
            except FileNotFoundError:
                self.console.print('============>File does not exist')
                continue
            except:
                self.console.print('============>Other errors...')
                continue
        self.pfwriter.generic(empty, self.sv_empty_fp + 'is_empty_table.txt')
        return 'Finished'


if __name__ == "__main__":
    from pypropel.path import to

    import pandas as pd

    prot_df = pd.DataFrame({
        'prot': ['1aig', '1aij', '1xqf'],
        'chain': ['L', 'L', 'A'],
    })

    p = IsEmpty(
        prot_df=prot_df,
        fasta_fp=to('data/fasta/'),
        sv_empty_fp=to('data/'),
    )
    print(p.fasta())