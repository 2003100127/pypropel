__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"


class Name:

    def __init__(self, ):
        pass

    def chain(self, prot_chain):
        return str(prot_chain) + 'l' if str(prot_chain).islower() else str(prot_chain)

    def seqchain(self, prot_chain):
        return str(prot_chain[0])