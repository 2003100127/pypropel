__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import os
import sys
sys.path.append(os.path.dirname(os.getcwd()) + '/')

import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from pypropel.prot.structure.distance import Distance
from pypropel.util.Console import Console


class AllAgainstAll(Distance.distance):

    def __init__(
            self,
            pdb_fp,
            pdb_name,
            verbose : bool = True,
    ):
        self.bio_parser = PDBParser()
        self.pdb_fpn = pdb_fp + pdb_name + '.pdb'
        self.structure = self.bio_parser.get_structure(pdb_name, self.pdb_fpn)
        self.model = self.structure[0]

        self.console = Console()
        self.console.verbose = verbose

    def chains(self, ):
        return [c.get_id() for c in self.structure.get_chains()]

    def calculate(self, ):
        chains = self.chains()
        num_chains = len(chains)
        df = pd.DataFrame()
        for i in range(num_chains):
            for j in range(i + 1, num_chains):
                self.console.print("============>chain 1 {} and chain2 {}".format(chains[i], chains[j]))
                working_chain1 = self.model[chains[i]]
                working_chain2 = self.model[chains[j]]
                # print(working_chain1, working_chain2)
                dist_mat = self.one2one_all(
                    working_chain1,
                    working_chain2
                )
                df_chain = pd.DataFrame(dist_mat)
                df_chain = df_chain.rename(columns={
                    0: 'res_fas_id1',
                    1: 'res1',
                    2: 'res_pdb_id1',
                    3: 'res_fas_id2',
                    4: 'res2',
                    5: 'res_pdb_id2',
                    6: 'dist',
                })
                df_chain['chain1'] = str(chains[i])
                df_chain['chain2'] = str(chains[j])
                df = pd.concat([df, df_chain], axis=0)
        return df


if __name__ == "__main__":
    from pypropel.path import to

    p = AllAgainstAll(
        pdb_fp=to('data/pdb/complex/'),
        pdb_name='1aij',
    )
    print(p.chains())
    # print(p.calculate())