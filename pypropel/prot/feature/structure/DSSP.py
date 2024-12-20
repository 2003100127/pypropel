__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import pandas as pd
from Bio.PDB import PDBParser as biopypdb
from Bio.PDB.DSSP import DSSP as biopythondssp
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Reader import Reader as pfreader
from pypropel.util.Console import Console


class DSSP:

    def __init__(
            self,
            prot_name,
            prot_chain,
            verbose: bool = True,
    ):
        self.prot_name = prot_name
        self.prot_chain = prot_chain
        self.pfreader = pfreader()
        self.pfwriter = pfwriter()
        self.console = Console()
        self.console.verbose = verbose

    def run(
            self,
            pdb_fp,
            sv_fp=None,
            mode='chain',
    ):
        """
        'A' means chain in the complex, 4 means the pdb position
        concerning the dssp key: ('A', (' ', 4, ' '))

        Parameters
        ----------
        prot_name
            prot name
        prot_chain
            prot chain
        mode
            chain or complex

        Returns
        -------
        """
        file_chain = chainname().chain(self.prot_chain)
        if mode == 'complex':
            pdb_fpn = pdb_fp + self.prot_name + '.pdb'
        else:
            pdb_fpn = pdb_fp + self.prot_name + file_chain + '.pdb'
        structure = biopypdb().get_structure(
            self.prot_name + self.prot_chain,
            pdb_fpn
        )
        model = structure[0]
        dssp = biopythondssp(
            model,
            pdb_fpn,
            dssp='mkdssp'
        )
        qualified = []
        pdb_ids = []
        fasta_ids = []
        rsa = []
        for v in list(dssp.keys()):
            if v[0] == self.prot_chain:
                qualified.append(v)
            else:
                continue
        for i, v in enumerate(qualified):
            rsa.append(dssp[v][3])
            fasta_ids.append(i+1)
            pdb_ids.append(v[1][1])
        # print(pdb_ids)
        # print(fasta_ids)
        # print(rsa)
        rsa_df = pd.DataFrame({
            'fasta_id': fasta_ids,
            'pdb_ids': pdb_ids,
            'rsa': rsa,
        })
        # print(rsa_df)
        if sv_fp != None:
            self.pfwriter.generic(
                df=rsa_df,
                sv_fpn=sv_fp + self.prot_name + file_chain + '.rsa',
                header=True
            )
        return rsa_df

    def access(
            self,
            rsa_fp
    ):
        file_chain = chainname().chain(self.prot_chain)
        df_rsa = self.pfreader.generic(
            rsa_fp + self.prot_name + file_chain + '.rsa',
            df_sep='\t',
            header=0,
        )
        df_rsa.columns = [
            'rsa_fas_id',
            'rsa_pdb_ids',
            'rsa_prob',
        ]
        df_rsa['rsa_prob'] = df_rsa['rsa_prob'].astype(float)
        return df_rsa


if __name__ == "__main__":
    from pypropel.path import to


    prot_df = pd.DataFrame({
        'prot': ['3pux', '3rko', '3udc', '3vr8', '4kjs', '4pi2', ],
        'chain': ['G', 'A', 'A', 'D', 'A', 'C', ],
    })
    for i in prot_df.index:
        print('No.{}: protein: {} chain: {}'.format(i + 1, prot_df.loc[i, 'prot'], prot_df.loc[i, 'chain']))
        p = DSSP(
            prot_name=prot_df.loc[i, 'prot'],
            prot_chain=prot_df.loc[i, 'chain'],
        )
        # p.run(
        #     pdb_fp='data/pdb/pdbtm/',
        #     sv_fp='data/rsa/',
        # )
        df_rsa = p.access(
            rsa_fp=to('data/rsa/')
        )
        print(df_rsa)