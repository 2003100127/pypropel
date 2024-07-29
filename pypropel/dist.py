__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List, Dict

import pandas as pd

from pypropel.prot.structure.distance.isite.heavy.AllAgainstAll import AllAgainstAll
from pypropel.prot.structure.distance.isite.heavy.OneToOne import OneToOne
from pypropel.prot.structure.distance.isite.heavy.Run import Run
from pypropel.prot.structure.distance.isite.check.Complex import Complex
from pypropel.prot.structure.distance.isite.check.Pair import Pair
from pypropel.prot.structure.distance.isite.Label import Label
from pypropel.prot.structure.distance.isite.check.TransmitterComplex import TransmitterComplex


def one_vs_one(
        pdb_path1 : str,
        pdb_name1 : str,
        file_chain1 : str,
        seq_chain1 : str,
        pdb_path2 : str,
        pdb_name2 : str,
        file_chain2 : str,
        seq_chain2 : str,
) -> Dict:
    return OneToOne(
        pdb_path1=pdb_path1,
        pdb_name1=pdb_name1,
        file_chain1=file_chain1,
        seq_chain1=seq_chain1,
        pdb_path2=pdb_path2,
        pdb_name2=pdb_name2,
        file_chain2=file_chain2,
        seq_chain2=seq_chain2,
    ).calculate()


def all_vs_all(
        pdb_fp,
        pdb_name,
):
    return AllAgainstAll(
        pdb_fp=pdb_fp,
        pdb_name=pdb_name,
    ).calculate()


def check_chain_complex(
        pdb_fp : str,
        prot_name : str,
        thres : float,
        sv_fp : str,
):
    return Complex(
        pdb_fp=pdb_fp,
        prot_name=prot_name,
        thres=thres,
        sv_fp=sv_fp,
    ).run()


def check_chain_paired(
        pdb_fp1 : str,
        pdb_fp2 : str,
        prot_name1 : str,
        prot_name2 : str,
        prot_chain1 : str,
        prot_chain2 : str,
        sv_fp : str,
        thres : float,
):
    return Pair(
        pdb_fp1=pdb_fp1,
        pdb_fp2=pdb_fp2,
        prot_name1=prot_name1,
        prot_name2=prot_name2,
        prot_chain1=prot_chain1,
        prot_chain2=prot_chain2,
        sv_fp=sv_fp,
        thres=thres,
    ).run()


def complex_calc_all(
        pdb_fp : str,
        prot_name : str,
        prot_chain : str,
        method : str,
        sv_fp : str,
):
    return Run(
        pdb_fp=pdb_fp,
        prot_name=prot_name,
        prot_chain=prot_chain,
        method=method,
        sv_fp=sv_fp,
    ).dist_without_aa()


def complex_calc_inter(
        pdb_fp : str,
        prot_name : str,
        prot_chain : str,
        method : str,
        sv_fp : str,
):
    return Run(
        pdb_fp=pdb_fp,
        prot_name=prot_name,
        prot_chain=prot_chain,
        method=method,
        sv_fp=sv_fp,
    ).dist_with_aa()


def cloud_check(
        order_list,
        job_fp,
        job_fn,
        cpu,
        memory,
        method,
        submission_method,
):
    return TransmitterComplex(
        order_list=order_list,
        job_fp=job_fp,
        job_fn=job_fn,
        cpu=cpu,
        memory=memory,
        method=method,
        submission_method=submission_method,
    ).execute()


def labelling(
        dist_fp,
        prot_name,
        file_chain,
        cutoff=6,
) -> pd.DataFrame:
    return Label(
        dist_fp=dist_fp,
        prot_name=prot_name,
        file_chain=file_chain,
        cutoff=cutoff,
    ).attach()


def interation_partners(
        dist_fp,
        prot_name,
        file_chain,
        pdb_fp,
        cutoff=6,
) -> List:
    return Label(
        dist_fp=dist_fp,
        prot_name=prot_name,
        file_chain=file_chain,
        cutoff=cutoff,
    ).partner(
        pdb_fp=pdb_fp,
    )


if __name__ == "__main__":
    from pypropel.path import to

    # dist_mat = one_vs_one(
    #     pdb_path1=to('data/pdb/complex/pdbtm/'),
    #     pdb_name1='1aij',
    #     file_chain1='',
    #     seq_chain1='L',
    #     pdb_path2=to('data/pdb/complex/pdbtm/'),
    #     pdb_name2='1aij',
    #     file_chain2='',
    #     seq_chain2='M',
    # )
    # df_dist = pd.DataFrame(dist_mat)
    # df_dist = df_dist.rename(columns={
    #     0: 'res_fas_id1',
    #     1: 'res1',
    #     2: 'res_pdb_id1',
    #     3: 'res_fas_id2',
    #     4: 'res2',
    #     5: 'res_pdb_id2',
    #     6: 'dist',
    # })
    # print(df_dist)

    # df_dist = all_vs_all(
    #     pdb_fp=to('data/pdb/complex/pdbtm/'),
    #     pdb_name='1aij',
    # )
    # print(df_dist)

    # print(check_chain_complex(
    #     pdb_fp=to('data/pdb/complex/pdbtm/'),
    #     prot_name='1aij',
    #     sv_fp=to('data/pdb/complex/pdbtm/'),
    #     thres=5.5,
    # ))

    # print(check_chain_paired(
    #     pdb_fp1=to('data/pdb/pdbtm/'),
    #     pdb_fp2=to('data/pdb/pdbtm/'),
    #     prot_name1='1aij',
    #     prot_name2='1aij',
    #     prot_chain1='L',
    #     prot_chain2='M',
    #     thres=6.,
    #     sv_fp=to('data/pdb/pdbtm/'),
    # ))

    # print(complex_calc_all(
    #     pdb_fp=to('data/pdb/complex/pdbtm/'),
    #     prot_name='1aij',
    #     prot_chain='L',
    #     method='heavy',
    #     sv_fp=to('data/pdb/complex/pdbtm/'),
    # ))

    # print(complex_calc_inter(
    #     pdb_fp=to('data/pdb/complex/pdbtm/'),
    #     prot_name='1aij',
    #     prot_chain='L',
    #     method='heavy',
    #     sv_fp=to('data/pdb/complex/pdbtm/'),
    # ))

    df_dist = labelling(
        dist_fp=to('data/pdb/complex/pdbtm/'),
        prot_name='1aij',
        file_chain='L',
        cutoff=6,
    )
    print(df_dist)

    print(interation_partners(
        dist_fp=to('data/pdb/complex/pdbtm/'),
        prot_name='1aij',
        file_chain='L',
        cutoff=6,
        pdb_fp=to('data/pdb/complex/pdbtm/'),
    ))
    ## ++++++++++++++++++++++++
    # from pypropel.util.Reader import Reader as pfreader
    #
    # df = pfreader().generic(df_fpn=to('data/ex/final.txt'), header=0)
    # prots = df.prot.unique()[2000:]
    #
    # param_config = {
    #     'pdb_fp': '-fp',
    #     'pdb_fn': '-fn',
    #     'sv_fp': '-op',
    # }
    # value_config = {
    #     'tool_fp': '/path/to/python',
    #     'script_fpn': './Complex.py',
    #     'pdb_fp': '/path/to/protein complex files/',
    #     'sv_fp': '/path/to/save/results/',
    # }
    #
    # for key, prot in enumerate(prots):
    #     order_list = [
    #         value_config['tool_fp'],
    #         value_config['script_fpn'],
    #
    #         param_config['pdb_fp'], value_config['pdb_fp'],
    #         param_config['pdb_fn'], prot,
    #         param_config['sv_fp'], value_config['sv_fp'],
    #     ]
    #     print(cloud_check(
    #         order_list=order_list,
    #         job_fp='/path/to/save/job files/',
    #         job_fn=str(key),
    #         cpu=2,
    #         memory=10,
    #         method='script',
    #         submission_method='sbatch',
    #     ))