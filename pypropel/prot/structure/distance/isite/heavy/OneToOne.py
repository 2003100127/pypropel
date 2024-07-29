__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from Bio.PDB.PDBParser import PDBParser
from pypropel.prot.structure.distance import Distance


class OneToOne(Distance.distance):

    def __init__(
            self,
            pdb_path1=None,
            pdb_name1=None,
            file_chain1=None,
            seq_chain1=None,
            pdb_path2=None,
            pdb_name2=None,
            file_chain2=None,
            seq_chain2=None,
    ):
        self.bio_parser = PDBParser()
        self.pdb_fpn1 = pdb_path1 + pdb_name1 + file_chain1 + '.pdb'
        self.pdb_fpn2 = pdb_path2 + pdb_name2 + file_chain2 + '.pdb'
        self.structure1 = self.bio_parser.get_structure(pdb_name1, self.pdb_fpn1)
        self.structure2 = self.bio_parser.get_structure(pdb_name2, self.pdb_fpn2)
        self.model1 = self.structure1[0]
        self.model2 = self.structure2[0]
        self.working_chain1 = self.model1[seq_chain1]
        self.working_chain2 = self.model2[seq_chain2]

    def calculate(self, ):
        return self.one2one_all(
            self.working_chain1,
            self.working_chain2,
        )


if __name__ == "__main__":
    from pypropel.path import to
    import pandas as pd

    p = OneToOne(
        pdb_path1=to('data/pdb/complex/pdbtm/'),
        pdb_name1='1aij',
        file_chain1='',
        seq_chain1='L',
        pdb_path2=to('data/pdb/complex/pdbtm/'),
        pdb_name2='1aij',
        file_chain2='',
        seq_chain2='M',
    )

    dist_mat = p.calculate()

    print(pd.DataFrame(dist_mat))