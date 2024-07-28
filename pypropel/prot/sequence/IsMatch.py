__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from pypropel.prot.sequence.Fasta import Fasta
from pypropel.prot.sequence.PDB import PDB
from pypropel.prot.sequence.XML import XML
from pypropel.util.Writer import Writer as pfwriter
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Console import Console


class IsMatch:

    def __init__(
            self,
            prot_df,
            fasta_path=None,
            pdb_path=None,
            xml_path=None,
            kind='pdb<->xml',
            sv_mismatch_fp='./',
            verbose: bool = True,
    ):
        self.prot_df = prot_df
        self.fasta_path = fasta_path
        self.pdb_path = pdb_path
        self.xml_path = xml_path
        self.kind = kind
        self.sv_mismatch_fp = sv_mismatch_fp

        self.pfwriter = pfwriter()
        self.console = Console()
        self.console.verbose = verbose

    def execute(self, ):
        mismatch = []
        for i, prot_name in enumerate(self.prot_df['prot']):
            prot_chain = self.prot_df['chain'][i]
            file_chain = chainname().chain(prot_chain)
            self.console.print('============>No{}. protein {} chain {}'.format(i, prot_name, prot_chain))
            try:
                if self.fasta_path:
                    fasta_seq = Fasta().get(
                        fasta_fpn=self.fasta_path + prot_name + file_chain + '.fasta'
                    )
                    # print(fasta_seq)
                else:
                    fasta_seq = None
                if self.pdb_path:
                    pdb_seq = PDB(
                        pdb_path=self.pdb_path,
                        pdb_name=prot_name,
                        file_chain=file_chain,
                        seq_chain=prot_chain,
                    ).chain()
                    # print(pdb_seq)
                else:
                    pdb_seq = None
                if self.xml_path:
                    xml_seq = XML().get(
                        xml_path=self.xml_path,
                        xml_name=prot_name,
                        seq_chain=prot_chain,
                    )
                    # print(xml_seq)
                else:
                    xml_seq = None
            except FileNotFoundError:
                self.console.print('============>File does not exist')
                continue
            except FileExistsError:
                self.console.print('============>File does not exist')
                continue
            except:
                self.console.print('============>Other errors...')
                continue
            else:
                if self.kind == 'fasta<->pdb':
                    if fasta_seq == pdb_seq:
                        self.console.print('============>They match.')
                    else:
                        self.console.print('============>They do not match.')
                        mismatch.append([prot_name, prot_chain])
                if self.kind == 'fasta<->xml':
                    if fasta_seq == xml_seq:
                        self.console.print('============>They match.')
                    else:
                        self.console.print('============>They do not match.')
                        mismatch.append([prot_name, prot_chain])
                if self.kind == 'pdb<->xml':
                    if pdb_seq == xml_seq:
                        self.console.print('============>They match.')
                    else:
                        self.console.print('============>They do not match.')
                        mismatch.append([prot_name, prot_chain])
        self.pfwriter.generic(mismatch, self.sv_mismatch_fp + 'is_match_table.txt')
        return 'Finished'


if __name__ == "__main__":
    from pypropel.path import to

    import pandas as pd
    prot_df = pd.DataFrame({
        'prot': ['1aig', '1aij', '1xqf'],
        'chain': ['L', 'L', 'A'],
    })

    p = IsMatch(
        prot_df=prot_df,
        fasta_path=to('data/fasta/'),
        pdb_path=to('data/pdb/pdbtm/'),
        xml_path=to('data/xml/'),
        kind='fasta<->xml',
        sv_mismatch_fp=to('data/'),
    )
    print(p.execute())