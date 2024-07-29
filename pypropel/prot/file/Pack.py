__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from pypropel.prot.structure.chain.Splitter import Splitter
from pypropel.prot.structure.chain.Format import Format
from pypropel.prot.structure.convert.ToFasta import ToFasta
from pypropel.prot.structure.hetatm.Remove import Remove as hetatmremover
from pypropel.prot.sequence.IsEmpty import IsEmpty
from pypropel.prot.sequence.IsMatch import IsMatch
from pypropel.util.FileIO import FileIO
from pypropel.util.Console import Console


class Pack:

    def __init__(
            self,
            prot_df,
            verbose: bool = True,
    ):
        self.prot_df = prot_df
        self.console = Console()
        self.console.verbose = verbose

    def execute(
            self,
            pdb_cplx_fp,
            pdb_fp,
            xml_fp,
            fasta_fp,
            kind='pdb<->xml',
    ):
        # ### /* block 1. split into chains */ ###
        self.console.print('=========>++++++++++++++++++++split into chains...\n+++++++++++++++++++++++++++')
        Splitter(
            prot_df=self.prot_df,
            pdb_path=pdb_cplx_fp,
            sv_fp=pdb_fp,
        ).pdb_per_chain()
        # ### /* block 2. delete END from PDB files */ ###
        self.console.print('=========>++++++++++++++++++++delete END from PDB files...\n+++++++++++++++++++++++++++')
        FileIO().makedir(pdb_fp + '/delend/')
        Format(
            prot_df=self.prot_df,
            sv_fp=pdb_fp + '/delend/',
        ).del_END_frompdb(
            pdb_path=pdb_fp,
        )
        # ### /* block 3. remove hetatm from PDB files */ ###
        self.console.print('=========>++++++++++++++++++++remove hetatm from PDB files...\n+++++++++++++++++++++++++++')
        hetatmremover(prot_df=self.prot_df).biopython(
            pdb_path=pdb_fp + '/delend/',
            sv_fp=pdb_fp,
        )
        # ### /* block 4. isMatch */ ###
        self.console.print('=========>++++++++++++++++++++is match...\n+++++++++++++++++++++++++++')
        IsMatch(
            prot_df=self.prot_df,
            pdb_path=pdb_fp + '/delend/',
            xml_path=xml_fp,
            sv_mismatch_fp=fasta_fp,
            kind=kind,
        ).execute()
        # ### /* block 5. ToFasta */ ###
        self.console.print('=========>++++++++++++++++++++to Fasta...\n+++++++++++++++++++++++++++')
        ToFasta(
            prot_df=self.prot_df,
            sv_fp=fasta_fp
        ).frompdb(
            pdb_path=pdb_fp + '/delend/',
        )
        # ### /* block 6. isEmpty */ ###
        self.console.print('=========>++++++++++++++++++++is empty...\n+++++++++++++++++++++++++++')
        IsEmpty(
            self.prot_df,
            sv_empty_fp=fasta_fp,
            fasta_fp=fasta_fp,
        ).fasta()
        return 'Finished'


if __name__ == "__main__":
    from pypropel.path import to

    import pandas as pd

    prot_df = pd.DataFrame({
        'prot': ['1aig', '1aij', '1xqf'],
        'chain': ['L', 'L', 'A'],
    })

    p = Pack(prot_df)

    print(p.execute(
        pdb_cplx_fp=to('data/pdb/complex/pdbtm/'),
        pdb_fp=to('data/'),
        xml_fp=to('data/xml/'),
        fasta_fp=to('data/'),
    ))