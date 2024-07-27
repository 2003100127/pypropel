__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from Bio import PDB


class Select(PDB.Select):

    def __init__(self, chain):
        self.chain_letters = chain

    def accept_chain(self, c):
        return (c.get_id() in self.chain_letters)

    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0
