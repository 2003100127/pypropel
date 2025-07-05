__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import os
from pypropel.util.FileIO import FileIO


class TransmitterComplex:

    def __init__(
            self,
            order_list,
            job_fp,
            job_fn,
            cpu,
            memory,
            method,
            submission_method,
    ):
        self.fileio = FileIO()
        self.order_list = order_list
        self.job_fp = job_fp
        self.job_fn = job_fn
        self.cpu = cpu
        self.memory = memory
        self.method = method
        self.submission_method = submission_method

    def initializer(self):
        switch = {
            "sbatch": self.sbatch,
            "qsub": self.qsub,
            "local": self.local,
        }
        execute_order = switch[self.submission_method]()
        print(f"======>The current order: {execute_order}")
        return execute_order

    def execute(self):
        """

        Returns
        -------

        """
        os.system(str(self.initializer()))
        return 'Submitted!'

    def sbatch(self, ):
        self.fileio.makedir(dir=self.job_fp)
        job_fpn = self.job_fp + self.job_fn + '.job'
        # print(job_fpn)
        with open(job_fpn, 'w') as job:
            # job.writelines("#!/bin/bash\n")
            job.writelines("#!/usr/bin/env bash\n")
            # job.writelines("#SBATCH --nodes=1\n")
            job.writelines("#SBATCH --time=23:00:00\n")
            job.writelines("#SBATCH --job-name={}\n".format(self.job_fn))
            job.writelines("#SBATCH --ntasks=1\n")
            job.writelines("#SBATCH --cpus-per-task={}\n".format(str(self.cpu)))
            job.writelines("#SBATCH --mem={}G\n".format(str(self.memory)))
            job.writelines("\n")
            job.writelines(str(self.construct_order())+"\n")
        return "sbatch {}".format(job_fpn)

    def qsub(self, ):
        pass

    def local(self, ):
        pass

    def construct_order(self, ):
        switch = {
            "script": self.script,
        }
        order = switch[self.method]()
        print(f"======>The constructed order is: {order}")
        return order

    def script(self, ):
        return self.recast(self.order_list)

    def recast(
            self,
            order_list,
            # send2cloud,
            # cloud_cmd,
    ):
        order = []
        for i in range(len(order_list)):
            order.append(order_list[i])
        a = " "
        if self.submission_method == 'qsub':
            # suffix = str(' "' + " | " + cloud_cmd)
            # b = a.join(order)
            # return "".join([b, suffix])
            return
        elif self.submission_method == 'sbatch':
            return a.join(order)
        else:
            return


if __name__ == "__main__":
    from pypropel.util.Reader import Reader as pfreader
    from pypropel.path import to

    # df = pfreader().generic(df_fpn='./cplx.txt')
    df = pfreader().generic(df_fpn='D:\Document\Programming\Python\minverse\minverse\data\deepisite/cplx.txt')
    # prots = df.prot.unique()[:2000]
    # prots = df.prot.unique()[2000:]
    prots = df[0].unique()
    print(prots)

    param_config = {
        'pdb_fp': '-fp',
        'pdb_fn': '-fn',
        'sv_fp': '-op',
        'thres': '-t',
    }
    value_config = {
        'tool_fp': '/project/cribbslab/jsun/miniconda3/envs/scflow/bin/python',
        'script_fpn': '/home/j/jsun/store/software/pypropel/pypropel/prot/structure/distance/isite/check/Complex.py',
        'pdb_fp': '/project/cribbslab/shared/proj048/sun/pdbtm/',
        'sv_fp': '/project/cribbslab/shared/proj048/sun/check/',
        'thres': '5.5',
    }

    for key, prot in enumerate(prots):
        order_list = [
            value_config['tool_fp'],
            value_config['script_fpn'],

            param_config['pdb_fp'], value_config['pdb_fp'],
            param_config['pdb_fn'], prot,
            param_config['thres'], value_config['thres'],
            param_config['sv_fp'], value_config['sv_fp'],
        ]
        p = TransmitterComplex(
            order_list=order_list,
            job_fp='/project/cribbslab/shared/proj048/sun/pdbtm/jobs/',
            job_fn=str(key),
            cpu=2,
            memory=10,
            method='script',
            submission_method='sbatch',
        )
        p.execute()