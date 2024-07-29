__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from pypropel.util.Reader import Reader as pfreader
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Console import Console


class Format:

    def __init__(
            self,
            prot_df,
            sv_fp,
            verbose: bool = True,
    ):
        self.prot_df = prot_df
        self.sv_fp = sv_fp

        self.pfreader = pfreader()
        self.console = Console()
        self.console.verbose = verbose

    def del_END_frompdb_op(
            self,
            pdb_path,
            prot_name,
            prot_chain,
    ):
        file_chain = chainname().chain(prot_chain)
        f = open(pdb_path + prot_name + file_chain + '.pdb')
        f_mirror = open(self.sv_fp + prot_name + file_chain + '.pdb', 'w')
        for line in f.readlines():
            # print(line)
            if not line.startswith('END'):
                f_mirror.write(line)
            else:
                continue
        # print(f)
        f_mirror.close()
        f.close()
        self.console.print('===============>Successfully reformatted')
        return 'Finished'

    def del_END_frompdb(
            self,
            pdb_path,
    ):
        for i in self.prot_df.index:
            prot_name = self.prot_df['prot'][i]
            prot_chain = self.prot_df['chain'][i]
            self.console.print('============>No{}. protein {} chain {}'.format(i, prot_name, prot_chain))
            self.del_END_frompdb_op(
                pdb_path=pdb_path,
                prot_name=prot_name,
                prot_chain=prot_chain,
            )
        return 'Finished'


if __name__ == "__main__":
    from pypropel.path import to

    import pandas as pd

    prot_df = pd.DataFrame({
        'prot': ['1aig', '1aij', '1xqf'],
        'chain': ['L', 'L', 'A'],
    })
    from pypropel.util.FileIO import FileIO

    p = Format(
        prot_df=prot_df,
        # sv_fp=FileIO().makedir(to('data/') + '/delend/'),
        sv_fp=to('data/'),
    )

    # print(p.del_END_frompdb_op(
    #     prot_name=DEFINE['prot_name'],
    #     prot_chain=DEFINE['prot_chain'],
    # ))

    print(p.del_END_frompdb(
        pdb_path=to('data/pdb/pdbtm/'),
    ))