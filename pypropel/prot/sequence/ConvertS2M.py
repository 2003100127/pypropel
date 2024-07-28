__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from pypropel.prot.sequence.Fasta import Fasta as sfasta
from pypropel.util.Reader import Reader as pfreader
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Console import Console


class ConvertS2M:

    def __init__(
            self,
            fasta_fp,
            prot_df,
            verbose: bool = True,
    ):
        self.fasta_fp = fasta_fp
        self.prot_df = prot_df

        self.sfasta = sfasta()
        self.cname = chainname()
        self.pfreader = pfreader()
        self.console = Console()
        self.console.verbose = verbose

    def integrate_seq(self, ):
        shape = self.prot_df.shape[0]
        rec = [[] for _ in range(shape)]
        self.console.print("=========>integrate {} protein sequences".format(shape))
        for i, prot_name in enumerate(self.prot_df['prot']):
            prot_chain = self.prot_df['chain'][i]
            file_chain = self.cname.chain(prot_chain)
            seq = self.sfasta.get(
                fasta_fpn=self.fasta_fp + prot_name + file_chain + '.fasta',
            )
            rec[i].append(prot_name + file_chain)
            rec[i].append(seq)
            # print(prot_name, prot_chain)
            # print(seq)
        return rec

    def save(
            self,
            list_2d,
            sv_fpn,
    ) -> str:
        self.console.print("=========>save fasta to files")
        len_list = len(list_2d)
        f = open(sv_fpn, 'w')
        for i in range(len_list):
            f.write('>' + str(list_2d[i][0])+'\n')
            f.write(str(list_2d[i][1])+'\n')
        f.close()
        self.console.print("=========>save finished")
        return 'Finished'


if __name__ == "__main__":
    from pypropel.path import to
    import pandas as pd

    p = ConvertS2M(
        fasta_fp=to('data/fasta/'),
        prot_df=pd.DataFrame({
            'prot': ['1aig', '1aij', '1xqf'],
            'chain': ['L', 'L', 'A'],
        })
    )

    seqs = p.integrate_seq()
    print(seqs)

    print(p.save(
        seqs,
        sv_fpn=to('data/fasta/s2m.fasta')
    ))