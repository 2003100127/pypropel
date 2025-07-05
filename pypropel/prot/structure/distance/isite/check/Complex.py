__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

# import os
# import sys
# sys.path.append(os.path.dirname(os.getcwd()) + '/')
import click
from pypropel.prot.structure.distance.isite.heavy.AllAgainstAll import AllAgainstAll as aaaheavy
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Console import Console


class Complex:

    def __init__(
            self,
            pdb_fp,
            prot_name,
            sv_fp,
            thres : float = 6,
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
        pypropel\pypropel>

        python ./prot/structure/distance/isite/check/Complex.py -fp ./data/pdb/complex/pdbtm/ -fn 1aij -op ./data/pdb/complex/pdbtm/

        Returns
        -------

        """
        self.console.print('=========>Protein PDB code: {}'.format(self.prot_name))
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
            self.console.print('=========>Protein chain 1: {}'.format(prot_chain1))
            chain1 = model[prot_chain1]
            for j in range(i+1, num_chains):
                prot_chain2 = chains[j]
                self.console.print('============>Protein chain 2: {}'.format(prot_chain2))
                chain2 = model[prot_chain2]
                if multimeric.check(chain1, chain2, thres=self.thres, verbose=self.verbose):
                    satisfied.append([self.prot_name, prot_chain1])
                    satisfied.append([self.prot_name, prot_chain2])
        return self.pfwriter.generic(
            satisfied,
            sv_fpn=self.sv_fp + self.prot_name + '.ccheck',
        )


@click.command(context_settings=dict(help_option_names=['-h', '--help']))
# @click.argument('method', type=str)
@click.option(
    '--pdb_fp', '-fp', type=str, required=True,
    help="""pdb file path"""
)
@click.option(
    '--pdb_fn', '-fn', type=str, required=True,
    help="""complex name"""
)
@click.option(
    '--thres', '-t', type=float, default=6.0, show_default=True,
    help="""threshold of dists"""
)
@click.option(
    '--sv_fp', '-o', type=str, default='./', show_default=True,
    help="""output path"""
)
@click.option(
    "--verbose", '-vb', default=True, show_default=True,
    help="whether to print output"
)
def cli(
    pdb_fp,
    pdb_fn,
    thres,
    sv_fp,
    verbose,
):
    """
    This function is used to check if two protein chains in the same protein complex are interacting.

    Examples
    --------
    pypropel_struct_check_cplx -fp D:/Document/Programming/Python/minverse/minverse/data/deepisite/pdbtm/cplx/ -fn 1aij -t 5.5 -o ./ -vb True

    """
    return Complex(
        pdb_fp=pdb_fp,
        prot_name=pdb_fn,
        thres=thres,
        sv_fp=sv_fp,
        verbose=verbose,
    ).run()


if __name__ == "__main__":
    cli()