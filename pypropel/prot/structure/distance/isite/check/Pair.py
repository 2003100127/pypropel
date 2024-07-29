__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import os
import sys

import pandas as pd

print(os.path.dirname(os.getcwd()) + '/')
sys.path.append(os.path.dirname(os.getcwd()) + '/')
from Bio.PDB.PDBParser import PDBParser
from pypropel.prot.structure.distance import Distance
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Console import Console


class Pair(Distance.distance):

    def __init__(
            self,
            pdb_fp1,
            pdb_fp2,
            prot_name1,
            prot_name2,
            prot_chain1,
            prot_chain2,
            sv_fp,
            thres : float = 6.,
            verbose : bool = True,
    ):
        self.pdb_fp1 = pdb_fp1
        self.pdb_fp2 = pdb_fp2
        self.prot_name1 = prot_name1
        self.prot_name2 = prot_name2
        self.prot_chain1 = prot_chain1
        self.prot_chain2 = prot_chain2
        self.sv_fp = sv_fp
        self.thres = thres
        self.file_chain1 = chainname().chain(self.prot_chain1)
        self.file_chain2 = chainname().chain(self.prot_chain2)
        self.bio_parser = PDBParser()
        self.pdb_fpn1 = self.pdb_fp1 + self.prot_name1 + self.file_chain1 + '.pdb'
        self.pdb_fpn2 = self.pdb_fp2 + self.prot_name2 + self.file_chain2 + '.pdb'
        self.structure1 = self.bio_parser.get_structure(self.prot_name1, self.pdb_fpn1)
        self.structure2 = self.bio_parser.get_structure(self.prot_name2, self.pdb_fpn2)
        self.model1 = self.structure1[0]
        self.model2 = self.structure2[0]
        self.working_chain1 = self.model1[self.prot_chain1]
        self.working_chain2 = self.model2[self.prot_chain2]

        self.pfwriter = pfwriter()
        self.console = Console()
        self.verbose = verbose
        self.console.verbose = self.verbose

    def calculate(self):
        pass

    def run(self, ):
        """

        Examples
        --------
        ./pypropel/pypropel>

        python ./prot/structure/distance/isite/check/Pair.py -fp1 ./data\pdb\pdbtm/ -fp2 ./data\pdb\pdbtm/ -fn1 1aij -fn2 1aij -c1 L -c2 M -t 6.0 -op ./data\pdb\pdbtm/

        Returns
        -------

        """
        self.console.print('=========>Protein PDB code 1: {}'.format(self.prot_name1))
        self.console.print('=========>Protein PDB chain 1: {}'.format(self.prot_chain1))
        self.console.print('=========>Protein PDB code 2: {}'.format(self.prot_name2))
        self.console.print('=========>Protein PDB chain 2: {}'.format(self.prot_chain2))
        ic = self.check(self.working_chain1, self.working_chain2, thres=self.thres, verbose=self.verbose)
        if ic:
            self.pfwriter.generic(
                df=pd.DataFrame([
                    [self.prot_name1, self.prot_chain1],
                    [self.prot_name2, self.prot_chain2],
                ]),
                sv_fpn=self.sv_fp + self.prot_name1 + self.prot_chain1 + '_' + self.prot_name2 + self.prot_chain2 + '.pcheck',
            )
        return ic


if __name__ == "__main__":
    # source = True
    source = False
    if source:
        import argparse
        parser = argparse.ArgumentParser(description='PPIs between a pair')

        parser.add_argument(
            "--pdb_fp1", "-fp1", dest='fp1', help='pdb file path 1', type=str
        )
        parser.add_argument(
            "--pdb_fp2", "-fp2", dest='fp2', help='pdb file path 2', type=str
        )
        parser.add_argument(
            "--pdb_fn1", "-fn1", dest='fn1', help='complex name 1', type=str
        )
        parser.add_argument(
            "--pdb_fn2", "-fn2", dest='fn2', help='complex name 2', type=str
        )
        parser.add_argument(
            "--pdb_c1", "-c1", dest='c1', help='complex chain 1', type=str
        )
        parser.add_argument(
            "--pdb_c2", "-c2", dest='c2', help='complex chain 2', type=str
        )
        parser.add_argument(
            "--thres", "-t", dest='t', help='threshold', type=float
        )
        parser.add_argument(
            "--sv_fp", "-op", dest='op', help='output path', type=str
        )
        args = parser.parse_args()
        if args.fp1:
            pdb_fp1 = args.fp1
        if args.fp2:
            pdb_fp2 = args.fp2
        if args.fn1:
            prot_name1 = args.fn1
        if args.fn2:
            prot_name2 = args.fn2
        if args.c1:
            prot_chain1 = args.c1
        if args.c2:
            prot_chain2 = args.c2
        if args.t:
            thres = args.t
        if args.op:
            sv_fp = args.op
    else:
        from pypropel.path import to

        pdb_fp1 = to('data/pdb/pdbtm/')
        pdb_fp2 = to('data/pdb/pdbtm/')
        prot_name1 = '1aij'
        prot_name2 = '1aij'
        prot_chain1 = 'L'
        prot_chain2 = 'M'
        thres = 6.
        sv_fp = to('data/pdb/pdbtm/')


    p = Pair(
        pdb_fp1=pdb_fp1,
        pdb_fp2=pdb_fp2,
        prot_name1=prot_name1,
        prot_name2=prot_name2,
        prot_chain1=prot_chain1,
        prot_chain2=prot_chain2,
        thres=thres,
        sv_fp=sv_fp,
    )

    p.run()