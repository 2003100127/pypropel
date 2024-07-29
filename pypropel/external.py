__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__developer__ = "Jianfeng Sun"
__maintainer__ = "Jianfeng Sun"
__email__="jianfeng.sunmt@gmail.com"

from typing import List, Dict

from pypropel.prot.feature.alignment.JSD import JSD


def jsd(
        order_list : List,
        job_fp,
        job_fn,
        cpu=2,
        memory=10,
        method='script',
        submission_method='local',
):
    return JSD(
        order_list=order_list,
        job_fp=job_fp,
        job_fn=job_fn,
        cpu=cpu,
        memory=memory,
        method=method,
        submission_method=submission_method,
    ).execute()


if __name__ == "__main__":
    from pypropel.path import to
    from pypropel.util.Reader import Reader as pfreader
    from pypropel.path import to

    # SR24_AtoI SR24_CtoU
    df = pfreader().generic(df_fpn=to('data/msa/clustal/wild/SR24_AtoI/prot.txt'))
    prots = df[0].unique()

    param_config = {
        'method': '-s',
        # 'window': '-w',
        # 'distance': '-d',
        'sv_fp': '-o',
        'clustal_fp': '',
    }
    value_config = {
        'tool_fp': 'python',
        'method': 'js_divergence',
        'window': '3',
        'distance': 'swissprot.distribution',
        'script_fpn': to('prot/feature/alignment/external/jsd/score_conservation.py'),
        'clustal_fp': to('data/msa/clustal/wild/SR24_AtoI/'),
        'sv_fp': to('data/jsd/SR24_AtoI/'),
    }

    for key, prot in enumerate(prots):
        order_list = [
            value_config['tool_fp'],
            value_config['script_fpn'],

            param_config['method'], value_config['method'],
            # param_config['window'], value_config['window'],
            # param_config['distance'], value_config['distance'],
            param_config['sv_fp'], value_config['sv_fp'] + prot + '.jsd',
            param_config['clustal_fp'], value_config['clustal_fp'] + prot + '.clustal',
        ]
        jsd(
            order_list=order_list,
            job_fp='./',
            job_fn=str(key),
        )
