__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import re
import pandas as pd
from Bio import SeqIO
from pypropel.util.Reader import Reader as pfreader
from pypropel.util.Console import Console


class ConvertM2S:

    def __init__(
            self,
            input_fpn : str ,
            in_format : str = 'fasta',
            sv_fp : str = './',
            verbose : bool = True,
    ):
        self.input_fpn = input_fpn
        self.in_format = in_format
        self.sv_fp = sv_fp

        self.fasta_ids = []
        self.fasta_seqs = []
        self.fasta_names = []
        self.fasta_dpts = []
        for i, record in enumerate(self.read()):
            self.fasta_ids.append(record.id)
            self.fasta_seqs.append(record.seq)
            self.fasta_names.append(record.name)
            self.fasta_dpts.append(record.description)
        self.df = pd.DataFrame({
            "fasta_ids": self.fasta_ids,
            "fasta_seqs": self.fasta_seqs,
            "fasta_names": self.fasta_names,
            "fasta_dpts": self.fasta_dpts,
        })
        self.df = self.drugbank_target_ids()
        self.df = self.drugbank_drug_ids()

        self.pfreader = pfreader()
        self.console = Console()
        self.console.verbose = verbose

        self.console.print("=========>content of the dataframe:\n {}".format(self.df.columns))
        self.console.print("=========>target IDs:\n {}".format(self.df.target_ids))

    def read(self, ):
        return SeqIO.parse(handle=self.input_fpn, format=self.in_format)

    def drugbank_target_ids(self, ):
        self.df['target_ids'] = self.df['fasta_ids'].apply(lambda x: re.sub(r'^.*\|', "", str(x)))
        # self.df['target_ids'] = self.df['fasta_ids'].apply(lambda x: re.split('\|', str(x))[1])
        return self.df

    def drugbank_drug_ids(self, ):
        self.df['drug_ids'] = self.df['fasta_dpts'].apply(lambda x: re.search(r'(?<=\()[^()]*(?=\)$)', str(x)).group().split('; '))
        return self.df

    def tofasta(self, ):
        self.console.print("=========>start to split into single fasta files".format(self.df.shape))
        self.console.print("=========># of proteins: {}".format(self.df.shape[0]))
        self.df = self.df.drop_duplicates(subset=['target_ids'], keep="first")
        self.console.print("=========># of proteins after deduplication: {}".format(self.df.shape[0]))
        for i in self.df.index:
            f = open(self.sv_fp + self.df.loc[i, 'target_ids'] + '.fasta', 'w')
            f.write('>' + self.df.loc[i, 'target_ids'] + '\n')
            f.write(str(self.df.loc[i, 'fasta_seqs']))
            f.close()
        return self.df


if __name__ == "__main__":
    from pypropel.path import to

    p = ConvertM2S(
        # input_fpn=to('data/msa/experimental_protein.fasta'),
        # in_format='fasta',
        # sv_fp=to('data/msa/'),
        input_fpn=to('data/yutong/protein.fasta'),
        in_format='fasta',
        sv_fp=to('data/yutong/'),
    )

    print(p.tofasta())
    print(p.df.target_ids)
    print(p.df.drug_ids)

    # from collections import Counter
    # print(Counter(p.drugbank_target_ids()))

    # print(p.drugbank_target_ids())
    # print(p.drugbank_drug_ids())

    # print(len(p.drugbank_target_ids()))
    # print(len(p.drugbank_drug_ids()))