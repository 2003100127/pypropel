__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np
import pandas as pd
from pypropel.prot.tool.isite.Reader import Reader as isitereader
from pypropel.util.Evaluate import Evaluate as mlgauge
from pypropel.util.Console import Console


class Segment():

    def __init__(
            self,
            prot_name,
            file_chain,
            dp,
            df_dist,
            pos_single_list,
            tool=None,
            tool_fp=None,
            sort_=None,
            verbose: bool = True,
    ):
        self.tool_fp = tool_fp
        self.prot_name = prot_name
        self.file_chain = file_chain
        self.pos_single_list = pos_single_list
        self.tool = tool
        self.sort_ = sort_
        self.mlevaluator = mlgauge()
        self.isitereader = isitereader()
        self.dp = dp
        self.df_dist = df_dist[[
            'fasta_id',
            'aa',
            'pdb_id',
            'is_contact'
        ]]
        self.row_real_dist = self.df_dist.shape[0]

        self.console = Console()
        self.console.verbose = verbose

    def fetch(self):
        switch = {
            'mbpred': self.mbpred,
            'delphi': self.delphi,
            'tma300': self.tma300,
        }
        return switch[self.tool]()

    def compare(self, target, cut_off, probability=False):
        # #/*** block 1. fetch target results ***/
        res_sorted = target.sort_values(['score'], ascending=False)
        # print(res_sorted)
        if probability:
            res_cutoff = res_sorted.loc[res_sorted['score'] > cut_off].reset_index(drop=True)
        else:
            res_cutoff = res_sorted.iloc[0: cut_off, :].reset_index(drop=True)
        # print(res_cutoff)
        res_cutoff = res_cutoff
        # print(res_cutoff)
        row_cutoff = res_cutoff.shape[0]
        row_offset = target.shape[0] - row_cutoff
        # #/*** block 1.1 y_true_cutoff ***/
        y_true_cutoff = []
        for i in range(row_cutoff):
            interid = self.df_dist['fasta_id'] == res_cutoff['fasta_id'][i]
            juery = (interid)
            tmp = self.df_dist.loc[juery]
            y_true_cutoff.append(tmp['is_contact'].sum())
        # print(y_true_cutoff)
        # #/*** block 1.2 y_true_all ***/
        y_true_all = list(res_sorted['is_contact'])
        # #/*** block 1.3 y_pred_cutoff ***/
        y_pred_cutoff = list((np.zeros([row_cutoff]) + 1).astype(np.int64))
        # print('y_pred: %s' % y_pred)
        # #/*** block 1.4 y_pred_all ***/
        y_pred_all = list(
            np.concatenate(
                [np.zeros([row_cutoff]) + 1, np.zeros([row_offset])],
                axis=0
            ).astype(np.int64)
        )
        # print(y_true_all, y_pred_all)
        # print(len(y_true_all), len(y_pred_all))
        # print('y_pred_all: %s' % y_pred_all)
        # #/*** block 1.5 y_score_all ***/
        y_score_all = list(res_sorted['score'])
        # print('y_score_all: %s' % y_score_all)

        # #/*** block 2. summary report ***/
        # pred_report = self.mlevaluator.all_metrics(y_true=venier, y_pred=y_pred)
        # pred_report = self.mlevaluator.all_metrics(y_true=real_con_all, y_pred=y_pred_all)
        # print(pred_report)

        # #/*** block 3. precision ***/
        precision = self.mlevaluator.precision(y_true=y_true_cutoff, y_pred=y_pred_cutoff)
        # print('------> precision: {}'.format(precision))
        # #/*** block 4. recall ***/
        recall = self.mlevaluator.recall(y_true=y_true_all, y_pred=y_pred_all)
        # print('------> recall: {}'.format(recall))
        # #/*** block 5. mcc ***/
        mcc = self.mlevaluator.mcc(y_true=y_true_all, y_pred=y_pred_all)
        # print('------> mcc: {}'.format(mcc))
        # #/*** block 6. f1score ***/
        f1score = self.mlevaluator.f1score(y_true=y_true_all, y_pred=y_pred_all)
        # print('------> f1score: {}'.format(f1score))
        # #/*** block 7. fbscore ***/
        fbscore = self.mlevaluator.fbscore(y_true=y_true_all, y_pred=y_pred_all, b=0.35)
#         # print('------> fbscore: {}'.format(fbscore))
        # #/*** block 8. accuracy ***/
        accuracy = self.mlevaluator.accuracy(y_true=y_true_all, y_pred=y_pred_all)
        # print('------> accuracy: {}'.format(accuracy))
        # #/*** block 9. hamming_loss ***/
        h_loss = self.mlevaluator.hamming_loss(y_true=y_true_all, y_pred=y_pred_all)
        # print('------> h_loss: {}'.format(h_loss))
        # #/*** block 10. balanced accuracy score ***/
        bal_acc = self.mlevaluator.balanced_accuracy_score(y_true=y_true_all, y_pred=y_pred_all)
        # print('------> bal_acc: {}'.format(bal_acc))
        # #/*** block 11. jaccard similarity score ***/
        jaccard = self.mlevaluator.jaccard_score(y_true=y_true_all, y_pred=y_pred_all)
        # print('------> jaccard: {}'.format(jaccard))
        # #/*** block 12. zero one loss ***/
        zo_loss = self.mlevaluator.zero_one_loss(y_true=y_true_all, y_pred=y_pred_all)
        # print('------> zo_loss: {}'.format(zo_loss))
        # #/*** block 13. confusion matrix ***/
        tp, fp, tn, fn = self.mlevaluator.cfmatrix1d(y_true=y_true_all, y_pred=y_pred_all)
#         # print('cf_matrix: {}'.format(tp, fp, tn, fn))
        # #/*** block 14. precision-recall ***/
        pr_precis, pr_recall, pr_thres = self.mlevaluator.prc(y_true=y_true_all, y_scores=y_score_all)
        ap = self.mlevaluator.ap(y_true=y_true_all, y_scores=y_score_all)
        # print('------> ap: {}'.format(ap))
        # #/*** block 15. roc ***/
        roc_fpr, roc_tpr, roc_thres = self.mlevaluator.roc(y_true=y_true_all, y_scores=y_score_all)
        auc = self.mlevaluator.auc(y_true=y_true_all, y_scores=y_score_all)
        # print('------> auc: {}'.format(auc))
        # #/*** block 16. precision-recall custom ***/
        pr_precis_custom, pr_recall_custom, pr_thres_custom = self.mlevaluator.prc_custom(y_true=y_true_all, y_scores=y_score_all)
        ap_custom = self.mlevaluator.ap(y_true=y_true_all, y_scores=y_score_all)
        # print('------> ap custom: {}'.format(ap_custom))
        # #/*** block 17. roc custom ***/
        roc_fpr_custom, roc_tpr_custom = self.mlevaluator.roc_custom(y_true=y_true_all, y_scores=y_score_all)
        auc_custom = self.mlevaluator.auc(y_true=y_true_all, y_scores=y_score_all)
        # print('------> auc custom: {}'.format(auc_custom))
        # #/*** block 17. roc custom ***/
        ppv = self.mlevaluator.ppv_custom(y_true=y_true_all, y_pred=y_pred_all)
        # print('------> ppv custom: {}'.format(ppv))
        # #/*** block 18. precision_by_hand ***/
        # precis_hand1 = sum(y_true_cutoff) / cut_off
        # precis_hand2 = res_cutoff['is_contact'].mean()
        # print('venier: %s' % venier)
        # print('------> precision1 man-made: {}'.format(precis_hand1))
        # print('------> precision2 man-made: {}'.format(precis_hand2))

        # #/*** block 18. plot ***/
        # plotsole().pr(pr_recall, pr_precis)
        # plotsole().roc(roc_fpr, roc_tpr)
        # plotsole().pr(pr_recall_custom, pr_precis_custom)
        # plotsole().roc(roc_fpr_custom, roc_tpr_custom)
        metrics_summary = {
            'precision': precision,
            'recall': recall,
            'f1score': f1score,
            'fbscore': fbscore,
            'mcc': mcc,
            'accuracy': accuracy,
            'h_loss': h_loss,
            'bal_acc': bal_acc,
            'jaccard': jaccard,
            'zo_loss': zo_loss,
            'tp': tp,
            'fp': fp,
            'tn': tn,
            'fn': fn,
            'pr_precis': pr_precis,
            'pr_recall': pr_recall,
            'pr_thres': pr_thres,
            'ap': ap,
            'roc_fpr': roc_fpr,
            'roc_tpr': roc_tpr,
            'roc_thres':roc_thres,
            'auc': auc,
            'pr_precis_custom': pr_precis_custom,
            'pr_recall_custom': pr_recall_custom,
            'roc_fpr_custom': roc_fpr_custom,
            'roc_tpr_custom': roc_tpr_custom,
            'ppv_custom': ppv,
            'y_true_all': y_true_all,
            'y_score_all': y_score_all
        }
        return metrics_summary

    def confer(self, y_true, y_pred):
        if y_true == 1 and y_pred == 1:
            return 'tp'
        elif y_true == 0 and y_pred == 1:
            return 'fp'
        elif y_true == 0 and y_pred == 0:
            return 'tn'
        else:
            return 'fn'

    def rsa(self, target, cut_off, rsa_path, probability=False, rsa_thres=[0.25]):
        res_sorted = target.sort_values(['score'], ascending=False)
        # print(res_sorted)
        if probability:
            res_cutoff = res_sorted.loc[res_sorted['score'] > cut_off].reset_index(drop=True)
        else:
            res_cutoff = res_sorted.iloc[0: cut_off, :].reset_index(drop=True)
        # print(res_cutoff)
        res_cutoff = res_cutoff
        # print(res_cutoff)
        row_cutoff = res_cutoff.shape[0]
        row_offset = target.shape[0] - row_cutoff
        y_true_cutoff = []
        for i in range(row_cutoff):
            interid = self.df_dist['fasta_id'] == res_cutoff['fasta_id'][i]
            juery = (interid)
            tmp = self.df_dist.loc[juery]
            y_true_cutoff.append(tmp['is_contact'].sum())
        # print(len(y_true_cutoff))
        y_pred_all = list(
            np.concatenate(
                [np.zeros([row_cutoff]) + 1, np.zeros([row_offset])],
                axis=0
            ).astype(np.int64)
        )
        res_sorted['y_pred'] = y_pred_all
        res_sorted['mark'] = res_sorted.apply(lambda x: self.confer(
            y_true=x['is_contact'],
            y_pred=x['y_pred'],
        ), axis=1)
        # print(res_sorted)
        res_gp = res_sorted.groupby(by=['mark'])
        res_gp_keys = res_gp.groups.keys()
        if 'fp' in res_gp_keys:
            res_gp_fp = res_gp.get_group('fp')
        else:
            res_gp_fp = pd.DataFrame(columns=res_sorted.columns)
        if 'fn' in res_gp_keys:
            res_gp_fn = res_gp.get_group('fn')
        else:
            res_gp_fn = pd.DataFrame(columns=res_sorted.columns)
        # print(res_gp_fp)
        # print(res_gp_fn)
        # print('y_score_all: %s' % y_score_all)
        rsa_df = self.isitereader.pfreader.generic(
            df_fpn=rsa_path + self.prot_name + self.file_chain + '.rsa',
            header=0,
        )
        # print(rsa_df)
        buried_fp_arr = []
        exposed_fp_arr = []
        buried_fn_arr = []
        exposed_fn_arr = []
        for rsa_thre in rsa_thres:
            rsa_df['state'] = rsa_df['rsa'].apply(
                lambda x: 'exposed' if x > rsa_thre else 'buried'
            )
            # print(rsa_df)
            fp_ids = res_gp_fp['fasta_id'].tolist()
            fn_ids = res_gp_fn['fasta_id'].tolist()
            # print(fp_ids)
            # print(fn_ids)
            buried_fp_num = rsa_df.loc[rsa_df['fasta_id'].isin(fp_ids)].loc[rsa_df['state'] == 'buried'].shape[0]
            exposed_fp_num = rsa_df.loc[rsa_df['fasta_id'].isin(fp_ids)].loc[rsa_df['state'] == 'exposed'].shape[0]
            buried_fn_num = rsa_df.loc[rsa_df['fasta_id'].isin(fn_ids)].loc[rsa_df['state'] == 'buried'].shape[0]
            exposed_fn_num = rsa_df.loc[rsa_df['fasta_id'].isin(fn_ids)].loc[rsa_df['state'] == 'exposed'].shape[0]
            # print(buried_fp_num)
            # print(exposed_fp_num)
            # print(buried_fn_num)
            # print(exposed_fn_num)
            buried_fp_arr.append(buried_fp_num)
            exposed_fp_arr.append(exposed_fp_num)
            buried_fn_arr.append(buried_fn_num)
            exposed_fn_arr.append(exposed_fn_num)
        buried_fp = ';'.join(np.array(buried_fp_arr).astype(np.str))
        exposed_fp = ';'.join(np.array(exposed_fp_arr).astype(np.str))
        buried_fn = ';'.join(np.array(buried_fn_arr).astype(np.str))
        exposed_fn = ';'.join(np.array(exposed_fn_arr).astype(np.str))
        buried = buried_fp + '+' + buried_fn
        exposed = exposed_fp + '+' + exposed_fn
        # print(buried)
        # print(exposed)
        # print(rsa_df)
        # tp, fp, tn, fn = self.mlevaluator.cfmatrix1d(y_true=y_true_all, y_pred=y_pred_all)
        # print('cf_matrix: {} {} {} {}'.format(tp, fp, tn, fn))
        metrics_summary = {
            'buried': buried,
            'exposed': exposed
        }
        return metrics_summary

    def mbpred(self):
        mbp_dist, chain_dist_all = self.isitereader.mbpred(
            mbp_path=self.tool_fp,
            file_name=self.prot_name,
            file_chain=self.file_chain,
            pos_single_list=self.pos_single_list,
            dp=self.dp,
            df_dist=self.df_dist,
            sort_=self.sort_
        )
        # chain_dist_all.col
        ensemble = pd.concat([mbp_dist, chain_dist_all], axis=1)
        # print(ensemble)
        # print(ensemble.shape[1])
        # print(ensemble.iloc[1])
        return ensemble

    def delphi(self):
        delphi_dist, chain_dist_all = self.isitereader.delphi(
            delphi_path=self.tool_fp,
            file_name=self.prot_name,
            file_chain=self.file_chain,
            pos_single_list=self.pos_single_list,
            dp=self.dp,
            df_dist=self.df_dist,
            sort_=self.sort_
        )
        # chain_dist_all.col
        ensemble = pd.concat([delphi_dist, chain_dist_all], axis=1)
        # print(ensemble)
        # print(ensemble.shape[1])
        # print(ensemble.iloc[1])
        return ensemble

    def tma300(self):
        tma300_dist, chain_dist_all = self.isitereader.tma300(
            tma300_path=self.tool_fp,
            file_name=self.prot_name,
            file_chain=self.file_chain,
            pos_single_list=self.pos_single_list,
            dp=self.dp,
            df_dist=self.df_dist,
            sort_=self.sort_
        )
        # print(tma300_dist)
        # print(chain_dist_all)
        ensemble = pd.concat([tma300_dist, chain_dist_all], axis=1)
        # print(ensemble)
        # print(ensemble.shape[1])
        # print(ensemble.iloc[1])
        return ensemble


if __name__ == "__main__":
    from pypropel.path import to
    from pypropel.prot.sequence.Fasta import Fasta as sfasta
    import tmkit as tmk
    from tmkit.seqnetrr.combo.Length import length as lscenario
    from tmkit.position.scenario.Segment import Segment as sscenario
    from pypropel.prot.structure.distance.isite.Label import Label as dlabel

    INIT = {
        'delphi_path': to('data/predictor/ppi/delphi/tm_alpha_n30/'),
        'mbp_cyto_path': to('data/predictor/ppi/mbpred/mbpredcyto/tm_alpha_n60/'),
        'mbp_tmh_path': to('data/predictor/ppi/mbpred/mbpredtm/tm_alpha_n60/'),
        'mbp_extra_path': to('data/predictor/ppi/mbpred/mbpredextra/tm_alpha_n60/'),
        # 'mbp_combined_path': to('data/predictor/ppi/mbpred/tm_alpha_n60/combined/phobius/'),
        'sv_pdbtm_fpn': to('data/predictor/ppi/mbpred/mbpredcombined/tm_alpha_n60/pdbtm/'),
        'sv_phobius_fpn': to('data/predictor/ppi/mbpred/mbpredcombined/tm_alpha_n60/phobius/'),
    }

    dp = dlabel(
        dist_fp=to('data/pdb/complex/pdbtm/'),
        prot_name='1aij',
        file_chain='L',
    )

    df_dist = dp.attach()
    print(df_dist)

    sequence = sfasta().get(
        fasta_fpn=to("data/fasta/1aijL.fasta")
    )
    print(sequence)
    len_seq = len(sequence)
    length_pos_list = lscenario().tosgl(len_seq)
    # print(length_pos_list)

    pdbtm_seg, pred_seg = tmk.topo.cepdbtm(
        pdb_fp=to('data/pdb/pdbtm/'),
        prot_name='1aij',
        seq_chain='L',
        file_chain='L',
        topo_fp=to('data/phobius/'),
        xml_fp=to('data/xml/'),
        fasta_fp=to('data/fasta/'),
    )
    print(pdbtm_seg)
    print(pred_seg)

    pdbtm_pos_cyto = sscenario().to_single(pdbtm_seg['cyto_lower'], pdbtm_seg['cyto_upper'])
    pdbtm_pos_tmh = sscenario().to_single(pdbtm_seg['tmh_lower'], pdbtm_seg['tmh_upper'])
    pdbtm_pos_extra = sscenario().to_single(pdbtm_seg['extra_lower'], pdbtm_seg['extra_upper'])
    phobius_pos_cyto = sscenario().to_single(pred_seg['cyto_lower'], pred_seg['cyto_upper'])
    phobius_pos_tmh = sscenario().to_single(pred_seg['tmh_lower'], pred_seg['tmh_upper'])
    phobius_pos_extra = sscenario().to_single(pred_seg['extra_lower'], pred_seg['extra_upper'])

    pdbtm_pos = pdbtm_pos_cyto + pdbtm_pos_tmh + pdbtm_pos_extra
    phobius_pos = phobius_pos_cyto + phobius_pos_tmh + phobius_pos_extra
    # ### /* ob init */
    print(phobius_pos)
    p = Segment(
        prot_name='1aij',
        file_chain='L',
        dp=dp,
        df_dist=df_dist,
        pos_single_list=phobius_pos,
        tool_fp=to('data/isite/deeptminter/'),
        # tool_fp=INIT['delphi_path'],
        # tool_fp=INIT['mbp_combined_pdbtm_path'],
        tool="tma300",
        sort_=1,
    )

    tool_results = p.fetch()
    # print(tool_results)

    comp = p.compare(
        target=tool_results,
        cut_off=147,
        probability=False

        # cut_off=0.5,
        # probability=True
    )
    # print(comp)

    # comp_rsa = p.rsa(
    #     target=tool_results,
    #     cut_off=5,
    #     probability=False,
    #     rsa_path=INIT['rsa_path'],
    #     rsa_thres=[0.1, 0.15, 0.2, 0.25],
    #     # cut_off=0.5,
    #     # probability=True
    # )
    # print(comp_rsa)