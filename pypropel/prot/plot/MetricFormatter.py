__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import re
import json
import functools
import numpy as np
import pandas as pd
from pypropel.util.Reader import Reader as pfreader
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Evaluate import Evaluate as mlgauge


class MetricFormatter:

    def __init__(
            self,
    ):
        self.pfreader = pfreader()
        self.pfwriter = pfwriter()

    def __get__(self, obj, type=None):
        pass

    def __call__(self, deal):
        @functools.wraps(deal)
        def switch(instance, *args, **kwargs):
            if kwargs['read_rest']:
                preds = self._read_rest(kwargs['pred_fpns'])
            if kwargs['read_cv']:
                preds['cv'] = self._read_cv(kwargs['cv_fpns'])
            return deal(instance, preds, **kwargs)
        return switch

    def _read_rest(self, pred_fpns, index=4):
        preds = {}
        for k, v in pred_fpns.items():
            pred = self.pfreader.generic(v).mean(axis=0).values.tolist()[:index]
            preds[k] = pred
        return preds

    def _read_cv(self, cv_fpns, index=4):
        preds = pd.DataFrame()
        for k, v in cv_fpns.items():
            preds[k] = self.pfreader.generic(v).mean(axis=0)
        return preds.mean(axis=1).values.tolist()[:index]

    def destring(
            self,
            roc_fpr,
    ):
        arr = [[] for _ in range(roc_fpr[0].shape[0])]
        for id, i in enumerate(roc_fpr[0]):
            i = re.sub(r'\]', '', i)
            i = re.sub(r'\[', '', i)
            h = i.split(', ')
            for j in h:
                arr[id].append(float(j))
        return arr

    def roc_custom_deprecate(self, fpr_fpn, tpr_fpn):
        roc_fpr = self.pfreader.generic(df_fpn=fpr_fpn)
        roc_tpr = self.pfreader.generic(df_fpn=tpr_fpn)
        # ###/*** block 1 ***/
        fpr = self.destring(roc_fpr)
        # print(fpr)
        # fpr_mean = list(np.mean(fpr, axis=0))
        fpr_mean = pd.DataFrame(fpr).mean(axis=0).tolist()
        # print(len(fpr_mean))
        fpr_mean = np.insert(fpr_mean, 0, 1)
        # ###/*** block 2 ***/
        tpr = self.destring(roc_tpr)
        # tpr_mean = np.mean(tpr, axis=0)
        tpr_mean = pd.DataFrame(tpr).mean(axis=0).tolist()
        tpr_mean = np.insert(tpr_mean, 0, 1)
        # plotsole().roc(fpr_mean, tpr_mean)
        return fpr_mean, tpr_mean

    def roc_json_one_prot(
            self,
            fpr_fpn : str,
            tpr_fpn : str,
    ):
        with open(fpr_fpn) as fp:
            fpr_dict = json.load(fp)
        with open(tpr_fpn) as fp:
            tpr_dict = json.load(fp)

        fpr_key_1st_outer = next(iter(fpr_dict))
        fpr_key_1st_inter = next(iter(fpr_dict[fpr_key_1st_outer]))
        # print(fpr_key_1st_outer)
        # print(fpr_key_1st_inter)
        fpr = fpr_dict[fpr_key_1st_outer][fpr_key_1st_inter]
        # print(fpr)
        tpr_key_1st_outer = next(iter(tpr_dict))
        tpr_key_1st_inter = next(iter(tpr_dict[tpr_key_1st_outer]))
        tpr = tpr_dict[tpr_key_1st_outer][tpr_key_1st_inter]
        fpr = np.insert(fpr, 0, 1).tolist()
        # print(fpr)
        tpr = np.insert(tpr, 0, 1).tolist()
        # print(tpr)
        # plotsole().roc(fpr_mean, tpr_mean)
        return fpr, tpr

    def pr_custom(self, p_fpn, r_fpn):
        p = self.pfreader.generic(df_fpn=p_fpn)
        r = self.pfreader.generic(df_fpn=r_fpn)
        precisions = self.destring(p)
        p_mean_ = list(np.mean(precisions, axis=0))
        # p_mean_.append(1)
        # p_mean_.insert(0, 0)
        recall = self.destring(r)
        r_mean_ = list(np.mean(recall, axis=0))
        # r_mean_.append(0)
        # r_mean_.insert(0, 1)
        # plotsole().roc(p_mean_, r_mean_)
        return p_mean_, r_mean_

    def roc_custom_all(self, y_true_fpn, y_score_fpn):
        y_trues = self.pfreader.generic(df_fpn=y_true_fpn).values.tolist()
        y_scores = self.pfreader.generic(df_fpn=y_score_fpn).values.tolist()
        fprs, tprs, thres = mlgauge().roc(y_true=y_trues, y_scores=y_scores)
        return fprs, tprs, thres

    def pr_custom_all(self, y_true_fpn, y_score_fpn):
        y_trues = self.pfreader.generic(df_fpn=y_true_fpn)
        y_scores = self.pfreader.generic(df_fpn=y_score_fpn)
        ps, rs, thres = mlgauge().prc(y_true=y_trues, y_scores=y_scores)
        return ps, rs, thres



if __name__ == "__main__":
    from pypropel.path import to

    p = MetricFormatter()
    # print(p.roc(
    #     fpr_fpn=to('data/al/prediction/rrc/tma165/tm_alpha_n57/tmh/pconsc4/pconsc4_roc_fpr.txt'),
    #     tpr_fpn=to('data/al/prediction/rrc/tma165/tm_alpha_n57/tmh/pconsc4/pconsc4_roc_tpr.txt'),
    # ))
    t1, t2 = p.roc_custom_deprecate(
        fpr_fpn=to('data/eval/tma300/tma300_roc_fpr_custom.txt'),
        tpr_fpn=to('data/eval/tma300/tma300_roc_tpr_custom.txt'),
    )
    print(t1)
    print(t2)

    t1, t2 = p.roc_json_one_prot(
        fpr_fpn=to('data/eval/tma300/tma300_roc_fpr_custom.json'),
        tpr_fpn=to('data/eval/tma300/tma300_roc_tpr_custom.json'),
    )
    print(t1)
    print(t2)

    # print(p.pr_custom(
    #     p_fpn=to('data/al/prediction/rrc/tma165/tm_alpha_n57/tmh/pconsc4/pconsc4_pr_precs.txt'),
    #     r_fpn=to('data/al/prediction/rrc/tma165/tm_alpha_n57/tmh/pconsc4/pconsc4_pr_recall.txt'),
    # ))