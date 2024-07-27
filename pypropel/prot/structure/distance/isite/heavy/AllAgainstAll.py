__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import os
import sys
sys.path.append(os.path.dirname(os.getcwd()) + '/')
from Bio.PDB.PDBParser import PDBParser
from pyprocpp.prot.structure.distance import Distance


class AllAgainstAll(Distance.distance):

    def __init__(
            self,
            pdb_fp,
            pdb_name,
    ):
        self.bio_parser = PDBParser()
        self.pdb_fpn = pdb_fp + pdb_name + '.pdb'
        self.structure = self.bio_parser.get_structure(pdb_name, self.pdb_fpn)
        self.model = self.structure[0]

    def chains(self, ):
        return [c.get_id() for c in self.structure.get_chains()]

    def calculate(self, ):
        chains = self.chains()
        num_chains = len(chains)
        for i in range(num_chains):
            for j in range(i + 1, num_chains):
                print(chains[i], chains[j])
                working_chain1 = self.model[chains[i]]
                working_chain2 = self.model[chains[j]]
                print(working_chain1, working_chain2)
                self.one2one_all(
                    working_chain1,
                    working_chain2
                )


if __name__ == "__main__":
    from pyprocpp.path import to

    p = AllAgainstAll(
        pdb_fp=to('data/pdb/complex/'),
        pdb_name='1aij',
    )
    print(p.chains())
    # print(p.calculate())