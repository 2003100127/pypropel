__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from Bio import PDB
from pypropel.prot.structure.residue.Select import select
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Console import Console


class Remove:

    def __init__(
            self, prot_df,
            verbose: bool = True,
    ):
        self.prot_df = prot_df
        self.io = PDB.PDBIO()
        self.parser = PDB.PDBParser()

        self.console = Console()
        self.console.verbose = verbose

    def biopython(
            self,
            pdb_path,
            sv_fp,
    ):
        """
        import subprocess

        for i, prot_name in enumerate(self.prot_df[0]):
            prot_chain = self.prot_df[1][i]
            if str(prot_chain).islower():
                c = str(prot_chain) + 'l'
            else:
                c = str(prot_chain)
            self.console.print('============>No.{}, {}'.format(i + 1, prot_name + prot_chain))
            order = str(
                'pdb_selchain -' +
                prot_chain.lower() + ' ' + pdb_path + prot_name + c +
                '.structure | pdb_delhetatm | pdb_tidy > ' +
                sv_fp + prot_name + c + '.structure'
            )
            self.console.print(order)
            s = subprocess.Popen(
                order,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE,
                shell=True
            )
            s.communicate()
        return 0
        Parameters
        ----------
        pdb_path
        sv_fp

        Returns
        -------

        """
        for i, prot_name in enumerate(self.prot_df['prot']):
            prot_chain = self.prot_df['chain'][i]
            file_chain = chainname().chain(prot_chain)
            self.console.print('============>No.{} protein {} chain {}'.format(i + 1, prot_name, prot_chain))
            try:
                struct = self.parser.get_structure(prot_name, pdb_path + prot_name + prot_chain + '.pdb')
                self.io.set_structure(struct)
                self.io.save(sv_fp + prot_name + file_chain + '.pdb', select=select())
            except:
                continue
        return 0


if __name__ == "__main__":
    from pypropel.path import to

    import pandas as pd

    prot_df = pd.DataFrame({
        'prot': ['1aig', '1aij', '1xqf'],
        'chain': ['L', 'L', 'A'],
    })

    p = Remove(prot_df=prot_df)

    print(p.biopython(
        pdb_path=to('data/pdb/'),
        sv_fp=to('data/pdb/'),
    ))