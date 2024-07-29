__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from pypropel.prot.sequence.PDB import PDB
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Console import Console


class ToFasta:
    
    def __init__(
            self,
            prot_df,
            sv_fp,
            verbose: bool = True,
    ):
        self.prot_df = prot_df
        self.sv_fp = sv_fp
        self.pfwriter = pfwriter()
        self.console = Console()
        self.console.verbose = verbose

    def frompdb(
            self,
            pdb_path,
    ):
        fails = []
        for i in self.prot_df.index:
            prot_name = self.prot_df.loc[i, 'prot']
            prot_chain = self.prot_df.loc[i, 'chain']
            file_chain = chainname().chain(prot_chain)
            self.console.print('============>No{}. protein {} chain {}'.format(i, prot_name, prot_chain))
            try:
                seq = PDB(
                    pdb_path=pdb_path,
                    pdb_name=prot_name,
                    file_chain=file_chain,
                    seq_chain=prot_chain,
                ).chain()
                if str(seq) == '':
                    continue
                else:
                    f = open(self.sv_fp + prot_name + file_chain + '.fasta', 'w')
                    f.write('>' + prot_name + prot_chain + '\n')
                    f.write(str(seq) + '\n')
                    f.close()
                    self.console.print('===============>successfully converted')
            except:
                fails.append([prot_name, prot_chain])
                print('No such a file: {}'.format(prot_name + prot_chain))
                continue
        self.pfwriter.generic(fails, self.sv_fp + 'fails_PDB2FASTA.txt')
        return 'finished'


if __name__ == "__main__":
    from pypropel.path import to

    import pandas as pd

    prot_df = pd.DataFrame({
        'prot': ['1aig', '1aij', '1xqf'],
        'chain': ['L', 'L', 'A'],
    })

    p = ToFasta(
        prot_df,
        sv_fp=to('data/'),
    )

    print(p.frompdb(
        pdb_path=to('data/pdb/pdbtm/'),
    ))