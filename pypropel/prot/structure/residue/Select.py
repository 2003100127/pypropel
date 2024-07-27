__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from Bio.PDB import Select


class select(Select):

    def __init__(self, ):
        pass

    def accept_residue(self, residue):
        return 1 if residue.id[0] == " " else 0