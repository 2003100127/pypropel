__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import mini3di
from Bio.PDB import PDBParser
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Console import Console


class Threedi:

    def __init__(
            self,
            pdb_fp,
            prot_name,
            prot_chain,
            verbose: bool = True,
    ):
        self.prot_name = prot_name
        self.prot_chain = prot_chain
        self.pdb_fp = pdb_fp
        self.threedi_encoder = mini3di.Encoder()

        self.console = Console()
        self.console.verbose = verbose

    def encode(
            self,
            mode='chain'
    ):
        """

        Parameters
        ----------
        mode
            chain or complex

        Returns
        -------

        """
        file_chain = chainname().chain(self.prot_chain)
        if mode == 'complex':
            pdb_fpn = self.pdb_fp + self.prot_name + '.pdb'
        else:
            pdb_fpn = self.pdb_fp + self.prot_name + file_chain + '.pdb'
        parser = PDBParser(QUIET=True)
        struct = parser.get_structure(self.prot_name + self.prot_chain, pdb_fpn)

        res = {}
        res[self.prot_name] = {}
        for chain in struct.get_chains():
            # print(chain.get_id())
            res[self.prot_name][chain.get_id()] = {}
            state = self.threedi_encoder.encode_chain(chain)
            encoded_sequence = self.threedi_encoder.build_sequence(state)
            # print(encoded_sequence)
            # print(state)
            res[self.prot_name][chain.get_id()]['state'] = state
            res[self.prot_name][chain.get_id()]['encoded_sequence'] = encoded_sequence
        # print(res)
        return res


if __name__ == "__main__":
    from pypropel.path import to
    import pandas as pd

    prot_df = pd.DataFrame({
        'prot': ['1aij', '1aig', '1xqf', ],
        'chain': ['L', 'L', 'A', ],
    })
    for i in prot_df.index:
        print('No.{}: protein: {} chain: {}'.format(i + 1, prot_df.loc[i, 'prot'], prot_df.loc[i, 'chain']))
        p = Threedi(
            prot_name=prot_df.loc[i, 'prot'],
            prot_chain=prot_df.loc[i, 'chain'],
            pdb_fp=to('data/pdb/pdbtm/'),
        )
        p.encode(mode='chain')

