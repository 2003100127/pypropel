__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import os
import sys
sys.path.append(os.path.dirname(os.getcwd()) + '/')
from pyprocpp.prot.structure.distance.isite.heavy.AllAgainstAll import AllAgainstAll as aaaheavy
from pyprocpp.util.Writer import Writer as pfwriter
from pyprocpp.util.Console import Console


class Complex:

    def __init__(
            self,
            pdb_fp,
            prot_name,
            sv_fp,
            thres=6,
            verbose: bool = True,
    ):
        self.pdb_fp = pdb_fp
        self.prot_name = prot_name
        self.sv_fp = sv_fp
        self.thres = thres

        self.pfwriter = pfwriter()
        self.console = Console()
        self.verbose = verbose
        self.console.verbose = self.verbose

    def run(self):
        """

        Examples
        --------
        pyprocpp\pyprocpp>

        python ./prot/structure/distance/isite/check/Complex.py -fp ./data/pdb/complex/pdbtm/ -fn 1aij -op ./data/pdb/complex/pdbtm/

        Returns
        -------

        """
        self.console.print('=========>Protein PDB code: {}'.format(prot_name))
        multimeric = aaaheavy(
            pdb_fp=self.pdb_fp,
            pdb_name=self.prot_name,
        )
        chains = multimeric.chains()
        num_chains = len(chains)
        model = multimeric.model
        satisfied = []
        for i in range(num_chains):
            prot_chain1 = chains[i]
            chain1 = model[prot_chain1]
            for j in range(i+1, num_chains):
                prot_chain2 = chains[j]
                chain2 = model[prot_chain2]
                if multimeric.check(chain1, chain2, thres=self.thres, verbose=self.verbose):
                    satisfied.append([self.prot_name, prot_chain1])
                    satisfied.append([self.prot_name, prot_chain2])
        return self.pfwriter.generic(
            satisfied,
            sv_fpn=self.sv_fp + self.prot_name + '.ccheck',
        )


if __name__ == "__main__":
    # source = True
    source = False
    if source:
        import argparse
        parser = argparse.ArgumentParser(description='PPIs in a complex')
        parser.add_argument(
            "--pdb_fp", "-fp", help='pdb file path', type=str
        )
        parser.add_argument(
            "--pdb_fn", "-fn", help='complex name', type=str
        )
        parser.add_argument(
            "--sv_fp", "-op", help='output path', type=str
        )
        args = parser.parse_args()
        # print(args)
        if args.pdb_fp:
            pdb_fp = args.pdb_fp
        if args.pdb_fn:
            prot_name = args.pdb_fn
        if args.sv_fp:
            sv_fp = args.sv_fp
    else:
        from pyprocpp.path import to

        pdb_fp = to('data/pdb/complex/pdbtm/')
        prot_name = '1aij'
        sv_fp = to('data/pdb/complex/pdbtm/')

    p = Complex(
        pdb_fp=pdb_fp,
        prot_name=prot_name,
        sv_fp=sv_fp,
    )
    p.run()