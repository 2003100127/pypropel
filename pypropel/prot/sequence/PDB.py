__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder


class PDB():

    def __init__(
            self,
            pdb_path,
            pdb_name,
            file_chain,
            seq_chain,
    ):
        self.pdb_path = pdb_path
        self.pdb_name = pdb_name
        self.file_chain = file_chain
        self.seq_chain = seq_chain
        self.pdb_fpn = self.pdb_path + self.pdb_name + self.file_chain + '.pdb'

    def chain(self, ):
        bio_parser = PDBParser()
        structure = bio_parser.get_structure(
            self.pdb_name,
            self.pdb_fpn,
        )
        model = structure[0]
        chain_init = model[self.seq_chain]
        ppb = PPBuilder()
        seq=[]
        for pp in ppb.build_peptides(chain_init):
            seq_tmp = str(pp.get_sequence())
            seq.append(seq_tmp)
        sequence = ''.join(seq)
        return sequence