__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from Bio import PDB
from pypropel.prot.structure.chain.Select import Select
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Console import Console


class Splitter:

    def __init__(
            self,
            prot_df,
            pdb_path=None,
            sv_fp=None,
            verbose: bool = True,
    ):
        self.prot_df = prot_df
        self.pdb_path = pdb_path
        self.sv_fp = sv_fp

        self.parser = PDB.PDBParser()
        self.io = PDB.PDBIO()

        self.pfwriter = pfwriter()
        self.chainname = chainname()
        self.console = Console()
        self.console.verbose = verbose

    def op(
            self,
            prot_name,
            prot_chain,
    ):
        file_chain = self.chainname.chain(prot_chain)
        seq_chain = self.chainname.seqchain(prot_chain)
        structure = self.parser.get_structure(
            prot_name,
            self.pdb_path + prot_name + '.pdb'
        )
        self.io.set_structure(structure)
        self.io.save(
            self.sv_fp + prot_name + file_chain + '.pdb',
            select=Select(seq_chain)
        )
        self.console.print('================>success in building ' + prot_name + prot_chain + ' model.')
        return 'Finished'

    def pdb_per_chain(
            self,
            mismatch_name='mismatch_records',
    ):
        """"""
        mismatch_records = []
        for (i, prot_name) in enumerate(self.prot_df['prot']):
            prot_chain = self.prot_df['chain'][i]
            self.console.print('============>No{}. protein {} chain {}'.format(i, prot_name, prot_chain))
            try:
                self.op(prot_name=prot_name, prot_chain=prot_chain)
            except:
                self.console.print('===============>mismatch_record: {}'.format([prot_name, prot_chain]))
                mismatch_records.append([prot_name, prot_chain])
                continue
        if self.sv_fp:
            self.pfwriter.generic(mismatch_records, sv_fpn=self.sv_fp + mismatch_name + '.txt')
        return 'Finished'

    def deprecated(
            self,
            pdb_path,
            sv_fp,
            sv_mismatch=True,
            mismatch_name='mismatch_records',
    ):
        """"""
        mismatch_records = []
        for i, prot_name in enumerate(self.prot_df[0]):
            print(i, prot_name)
            prot_chain = self.prot_df[1][i]
            try:
                pdb = self.parser.get_structure(prot_name, pdb_path + prot_name + ".pdb")
                count = 0
                for chain in pdb.get_chains():
                    # print(chain)
                    if chain.id == prot_chain:
                        self.io.set_structure(chain)
                        if str(prot_chain).islower():
                            file_chain = str(prot_chain) + 'l'
                        else:
                            file_chain = str(prot_chain)
                        self.io.save(sv_fp + pdb.get_id() + file_chain + ".pdb")
                        count = count + 1
                    else:
                        continue
                if count == 0:
                    print('mismatch_record: {}' .format([prot_name, prot_chain]))
                    mismatch_records.append([prot_name, prot_chain])
            except:
                continue
        if sv_mismatch:
            self.pfwriter.generic(mismatch_records, sv_fpn=sv_fp + mismatch_name)
        return True


if __name__ == "__main__":
    from pypropel.path import to

    import pandas as pd

    prot_df = pd.DataFrame({
        'prot': ['1aig', '1aij', '1xqf'],
        'chain': ['L', 'L', 'A'],
    })

    p = Splitter(
        prot_df=prot_df,
        pdb_path=to('data/pdb/complex/pdbtm/'),
        sv_fp=to('data/'),
    )

    print(p.pdb_per_chain())