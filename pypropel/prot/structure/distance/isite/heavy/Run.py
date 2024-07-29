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
from pypropel.prot.structure.distance.isite.heavy.AllAgainstAll import AllAgainstAll as aaaheavy
from pypropel.util.Writer import Writer as pfwriter
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Console import Console


class Run:

    def __init__(
            self,
            pdb_fp,
            prot_name,
            prot_chain,
            method,
            sv_fp,
            verbose: bool = True,
    ):
        self.pdb_fp = pdb_fp
        self.prot_name = prot_name
        self.prot_chain = prot_chain
        self.method = method
        self.sv_fp = sv_fp

        self.pfwriter = pfwriter()
        self.console = Console()
        self.verbose = verbose
        self.console.verbose = self.verbose

    def dist_without_aa(self, ):
        """
        order = str(
                'echo ' + '"' + 'conda activate common && python ' +
                self.cloud_fp + ' -pn ' + prot_name + ' -pc1 ' +
                prot_chain1 + ' -pdbp ' + self.pdb_path +
                ' -k heavy' + ' -sp ' + self.sv_fp + ' "' +
                " | qsub -l vf=2g -pe serial 4 -q all.q -N 'j.sun' "
            )
        Returns
        -------

        """
        self.console.print('=========>Protein PDB code: {}'.format(self.prot_name))
        self.console.print('=========>Chain of focus: {}'.format(self.prot_chain))
        multimeric = aaaheavy(
            pdb_fp=self.pdb_fp,
            pdb_name=self.prot_name,
        )
        chains = multimeric.chains()
        model = multimeric.model
        chains.remove(self.prot_chain)
        df_dist = pd.DataFrame()
        chain = model[self.prot_chain]
        pos_list = []
        for j, prot_chain2 in enumerate(chains):
            if prot_chain2 != self.prot_chain:
                self.console.print('============>Partner chain: {}'.format(prot_chain2))
                chain2 = model[chains[j]]
                if self.method == 'heavy':
                    dists = multimeric.one2one_minimal(chain, chain2, verbose=self.verbose)
                    df_dist[prot_chain2] = pd.DataFrame(dists)[3]
                    if j == 0:
                        for pos in dists:
                            pos_list.append(pos[:3])
                elif self.method == 'comprehensive':
                    pass
                else:
                    pass
        df = pd.concat([pd.DataFrame(pos_list), df_dist], axis=1)
        print(df)
        file_chain = chainname().chain(self.prot_chain)
        self.pfwriter.generic(
            df,
            self.sv_fp + self.prot_name + file_chain + '.dist',
            header=True,
        )
        return 0

    def dist_with_aa(self, thres=6):
        """
        order = str(
                'echo ' + '"' + 'source activate common && python ' +
                self.cloud_fp + ' -pn ' + prot_name + ' -pc1 ' +
                prot_chain1 + ' -pdbp ' + self.pdb_path +
                ' -k heavy' + ' -sp ' + self.sv_fp + ' "' +
                " | qsub -l vf=2g -pe serial 4 -q all.q -N 'j.sun' "
            )
        ..  @example:
            ---------
            python ./template/Mutual.py -pn 1fft -pc1 A -pdbp ./ -k heavy -sp ./
        :return:
        """
        self.console.print('=========>Protein PDB code: {}'.format(self.prot_name))
        self.console.print('=========>Chain of focus: {}'.format(self.prot_chain))
        file_chain1 = chainname().chain(self.prot_chain)
        multimeric = aaaheavy(
            pdb_fp=self.pdb_fp,
            pdb_name=self.prot_name,
        )
        chains = multimeric.chains()
        model = multimeric.model
        chains.remove(self.prot_chain)
        chain1 = model[self.prot_chain]
        accumulate = []
        for j, prot_chain2 in enumerate(chains):
            self.console.print('============>Partner chain: {}'.format(prot_chain2))
            chain2 = model[chains[j]]
            if self.method == 'heavy':
                dists = multimeric.one2one_all(chain1, chain2, verbose=self.verbose)
                df_dist = pd.DataFrame(dists)
                df_dist = df_dist.loc[df_dist[6] < thres].reset_index(drop=True)
                df_dist[7] = self.prot_chain
                df_dist[8] = prot_chain2
                df_dist = df_dist[[7, 2, 1, 8, 5, 4, 6]]
                df_dist.columns = [
                    'chain1',
                    'pos1',
                    'aa1',
                    'chain2',
                    'pos2',
                    'aa2',
                    'dist',
                ]
                accumulate.append(df_dist)
            elif self.method == 'comprehensive':
                pass
            else:
                pass
        df = pd.concat(accumulate, axis=0)
        self.pfwriter.generic(
            df=df,
            sv_fpn=self.sv_fp + self.prot_name + file_chain1 + '.pdbinter',
        )
        return


if __name__ == "__main__":
    # source = True
    source = False
    if source:
        import argparse
        parser = argparse.ArgumentParser(description='Distance of interaction sites')
        parser.add_argument(
            "--pdb_fp", "-fp", help='pdb file path', type=str
        )
        parser.add_argument(
            "--pdb_fn", "-fn", help='complex name', type=str
        )
        parser.add_argument(
            "--chain", "-c", help='protein chain of focus', type=str
        )
        parser.add_argument(
            "--method", "-m", help='method', type=str
        )
        parser.add_argument(
            "--sv_fp", "-op", help='output path', type=str
        )
        args = parser.parse_args()
        if args.pdb_fp:
            pdb_fp = args.pdb_fp
        if args.pdb_fn:
            prot_name = args.pdb_fn
        if args.chain:
            prot_chain = args.chain
        if args.method:
            method = args.method
        if args.sv_fp:
            sv_fp = args.sv_fp
    else:
        from pypropel.path import to

        pdb_fp = to('data/pdb/complex/pdbtm/')
        prot_name = '1aij'
        prot_chain = 'L'
        method = 'heavy'
        sv_fp = to('data/pdb/complex/pdbtm/')

    p = Run(
        pdb_fp=pdb_fp,
        prot_name=prot_name,
        prot_chain=prot_chain,
        method=method,
        sv_fp=sv_fp,
    )

    p.dist_without_aa()

    # p.dist_with_aa()