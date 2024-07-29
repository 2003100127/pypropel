__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import time
import pandas as pd
from pypropel.util.ComputLib import ComputLib as computlib
from pypropel.prot.structure.distance.isite.heavy.AllAgainstAll import AllAgainstAll as aaaheavy
from pypropel.util.Reader import Reader as pfreader
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Console import Console


class Label:

    def __init__(
            self,
            dist_fp,
            prot_name,
            file_chain,
            cutoff=6,
            verbose: bool = True,
    ):
        self.prot_name = prot_name
        self.file_chain = file_chain
        self.dist_fpn = dist_fp + self.prot_name + self.file_chain + '.dist'
        self.cutoff = cutoff

        self.pfreader = pfreader()
        self.computlib = computlib()
        self.console = Console()
        self.console.verbose = verbose

    def attach(self, ):
        self.console.print('================>Labeling data...')
        start_time = time.time()
        df_dist = self.pfreader.generic(self.dist_fpn, header=0)
        # print(df_dist)
        dists = df_dist.iloc[:, 3:]
        dist_mins = dists.min(axis=1)
        # print(dist_mins)
        inter_ids = dist_mins.loc[dist_mins < self.cutoff].index.tolist()
        noninter_ids = dist_mins.loc[dist_mins >= self.cutoff].index.tolist()
        df_dist['is_contact'] = -1
        df_dist.loc[inter_ids, 'is_contact'] = 1
        df_dist.loc[noninter_ids, 'is_contact'] = 0
        columns = ['fasta_id', 'aa', 'pdb_id']
        for i in range(dists.shape[1]):
            columns.append('dist_' + str(i+1))
        columns.append('is_contact')
        df_dist.columns = columns
        end_time = time.time()
        self.console.print('================>Time to label distances {} {}: {}s.'.format(self.prot_name, self.file_chain, end_time - start_time))
        return df_dist

    def segment(
            self,
            dist_df,
            pos_df,
            by : str | float = 'interact_id',
    ):
        dist_df_ = dist_df[[
            'fasta_id',
            'aa',
            'pdb_id',
            'is_contact',
        ]]
        num_samples = pos_df.shape[0]
        dists__ = pd.DataFrame()
        dists__[0] = dist_df_['fasta_id']
        dists__[1] = dist_df_.index.values
        dists__ = dists__.values.tolist()
        dist_dict = self.computlib.tactic6(dists__)
        dist_ids = []
        for i in range(num_samples):
            id = pos_df[by][i]
            # print(id)
            dist_id = dist_dict[id]
            dist_ids.append(dist_id)
        dist_df_ = dist_df_.iloc[dist_ids]
        dist_df_.columns = [
            'fasta_id',
            'aa',
            'pdb_id',
            'is_contact',
        ]
        dist_df_ = dist_df_.reset_index(inplace=False, drop=True)
        return dist_df_

    def partner(self, pdb_fp):
        multimeric = aaaheavy(
            pdb_fp=pdb_fp,
            pdb_name=self.prot_name,
        )
        chains = multimeric.chains()
        seq_chain = chainname().seqchain(self.file_chain)
        chains.remove(seq_chain)
        dist_df = self.pfreader.generic(self.dist_fpn, header=0)
        dists = dist_df.iloc[:, 3:]
        partners = []
        for i, col in enumerate(dists.columns):
            if dists[col].min() < self.cutoff:
                partners.append(chains[i])
        return partners


if __name__ == "__main__":
    from pypropel.path import to

    p = Label(
        dist_fp=to('data/pdb/complex/pdbtm/'),
        prot_name='1aij',
        file_chain='L',
        cutoff=6,
    )
    dists = p.attach()
    # print(dists)
    # print(dists.to_csv(to('data/pdb/complex/pdbtm/dist_label.txt')))

    # partners = p.partner(
    #     pdb_fp=to('data/pdb/complex/pdbtm/'),
    # )
    # print(partners)

    import tmkit as tmk

    lower_ids, upper_ids = tmk.topo.from_pdbtm(
        xml_fp=to('data/xml/'),
        prot_name='1aij',
        seq_chain='L',
        topo='tmh',
    )
    print(lower_ids)
    print(upper_ids)

    pdb_ids = tmk.seq.pdbid(
        pdb_fp=to('data/pdb/pdbtm/'),
        prot_name='1aij',
        seq_chain='L',
        file_chain='L',
    )
    print(pdb_ids)

    fas_ids = tmk.seq.fasid(
        fasta_fpn=to('data/fasta/1aijL.fasta'),
    )
    print(fas_ids)

    from tmkit.topology.pdbtm.ToFastaId import toFastaId

    fasta_lower_tmh, fasta_upper_tmh = toFastaId().tmh(
        pdbid_map=pdb_ids,
        fasid_map=fas_ids,
        xml_fp=to('data/xml/'),
        prot_name='1aij',
        seq_chain='L',
    )
    print(fasta_lower_tmh)
    print(fasta_upper_tmh)
    from tmkit.position.scenario.Segment import Segment as ppssegment

    pos_list_single = ppssegment().to_single(fasta_lower_tmh, fasta_upper_tmh)
    print(pos_list_single)

    df_seg = p.segment(
        dist_df=dists,
        pos_df=pd.DataFrame(pos_list_single),
        by=0,
    )
    print(df_seg)
    # print(df_seg.to_csv('./asdasd.txt'))