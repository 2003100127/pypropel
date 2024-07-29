__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import json
import pandas as pd
from pypropel.prot.tool.isite.Segment import Segment as sprecision
from pypropel.prot.structure.distance.isite.Label import Label as dlabel
from pypropel.prot.feature.sequence.Length import Length as fealength
from pypropel.prot.sequence.Name import Name as chainname
from pypropel.util.Reader import Reader as pfreader
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Console import Console


class Dispatcher:

    def __init__(
            self,
            verbose: bool = True,
    ):
        self.pfreader = pfreader()
        self.pfwriter = pfwriter()
        self.console = Console()
        self.console.verbose = verbose

    def segment(
            self,
            prot_df,
            dist_fp,
            dist_limit,
            tool_fp,
            tool,
            pdb_fp,
            topo_fp,
            xml_fp,
            fasta_fp,
            segment,
            sort,
            sv_fp=None,
    ):
        cutoffs = [1, 2, 5, 10]

        precis = {}
        recall = {}
        f1 = {}
        fb = {}
        mcc = {}
        acc = {}
        h_loss = {}
        bal_acc = {}
        jaccard = {}
        zo_loss = {}
        tp = {}
        fp = {}
        tn = {}
        fn = {}
        ppv_custom = {}
        auc = {}
        ap = {}

        roc_fpr = {}
        roc_tpr = {}
        roc_thres = {}
        roc_fpr_custom = {}
        roc_tpr_custom = {}

        pr_precs = {}
        pr_recall = {}
        pr_thres = {}
        pr_precs_custom = {}
        pr_recall_custom = {}

        y_true_all = []
        y_score_all = []
        for i, prot_name in enumerate(prot_df['prot']):
            prot_chain = prot_df.loc[i, 'chain']
            file_chain = chainname().chain(prot_chain)
            # print('No.%d: target: %s %s' % (i + 1, prot_name, file_chain))
            # /* distance dataframe */
            dp = dlabel(
                dist_fp=dist_fp,
                prot_name=prot_name,
                file_chain=file_chain,
                cutoff=dist_limit,
                verbose=True,
            )
            df_dist = dp.attach()
            pos_single_list, len_seg = self.pos_sgl_list(
                prot_name=prot_name,
                seq_chain=prot_chain,
                pdb_fp=pdb_fp,
                topo_fp=topo_fp,
                xml_fp=xml_fp,
                fasta_fp=fasta_fp,
                segment=segment,
            )
            precis_cutoff = {}
            recall_cutoff = {}
            f1_cutoff = {}
            fb_cutoff = {}
            mcc_cutoff = {}
            acc_cutoff = {}
            h_loss_cutoff = {}
            bal_acc_cutoff = {}
            jaccard_cutoff = {}
            zo_loss_cutoff = {}
            tp_cutoff = {}
            fp_cutoff = {}
            tn_cutoff = {}
            fn_cutoff = {}
            ppv_cutoff = {}
            auc_cutoff = {}
            ap_cutoff = {}

            roc_fpr_cutoff = {}
            roc_tpr_cutoff = {}
            roc_thres_cutoff = {}
            roc_fpr_custom_cutoff = {}
            roc_tpr_custom_cutoff = {}

            pr_precis_cutoff = {}
            pr_recall_cutoff = {}
            pr_thres_cutoff = {}
            pr_precis_custom_cutoff = {}
            pr_recall_custom_cutoff = {}

            y_true_all_cutoff = []
            y_score_all_cutoff = []
            # ###/*** block ***/
            if pos_single_list == []:
                continue
            else:
                init_compare = sprecision(
                    prot_name=prot_name,
                    file_chain=file_chain,
                    dp=dp,
                    df_dist=df_dist,
                    pos_single_list=pos_single_list,
                    tool_fp=tool_fp,
                    tool=tool,
                    sort_=sort
                )
                res_tool = init_compare.fetch()
                for j in range(len(cutoffs)):
                    cutoff = round(len_seg / cutoffs[j])
                    # print('---> Cutoff value is %d ***' % (cutoff))
                    metrics_summary = init_compare.compare(res_tool, cutoff)
                    precis_cutoff[cutoff] = metrics_summary['precision']
                    recall_cutoff[cutoff] = metrics_summary['recall']
                    f1_cutoff[cutoff] = metrics_summary['f1score']
                    fb_cutoff[cutoff] = metrics_summary['fbscore']
                    mcc_cutoff[cutoff] = metrics_summary['mcc']
                    acc_cutoff[cutoff] = metrics_summary['accuracy']
                    h_loss_cutoff[cutoff] = metrics_summary['h_loss']
                    bal_acc_cutoff[cutoff] = metrics_summary['bal_acc']
                    jaccard_cutoff[cutoff] = metrics_summary['jaccard']
                    zo_loss_cutoff[cutoff] = metrics_summary['zo_loss']
                    tp_cutoff[cutoff] = int(metrics_summary['tp'])
                    fp_cutoff[cutoff] = int(metrics_summary['fp'])
                    tn_cutoff[cutoff] = int(metrics_summary['tn'])
                    fn_cutoff[cutoff] = int(metrics_summary['fn'])
                    ppv_cutoff[cutoff] = metrics_summary['ppv_custom']

                    auc_cutoff[cutoff] = metrics_summary['auc']
                    ap_cutoff[cutoff] = metrics_summary['ap']

                    roc_fpr_cutoff[cutoff] = metrics_summary['roc_fpr'].tolist()
                    roc_tpr_cutoff[cutoff] = metrics_summary['roc_tpr'].tolist()
                    roc_thres_cutoff[cutoff] = metrics_summary['roc_thres'].tolist()
                    roc_fpr_custom_cutoff[cutoff] = metrics_summary['roc_fpr_custom']
                    roc_tpr_custom_cutoff[cutoff] = metrics_summary['roc_tpr_custom']

                    pr_precis_cutoff[cutoff] = metrics_summary['pr_precis'].tolist()
                    pr_recall_cutoff[cutoff] = metrics_summary['pr_recall'].tolist()
                    pr_thres_cutoff[cutoff] = metrics_summary['pr_thres'].tolist()
                    pr_precis_custom_cutoff[cutoff] = metrics_summary['pr_precis_custom']
                    pr_recall_custom_cutoff[cutoff] = metrics_summary['pr_recall_custom']
                    if j == 0:
                        y_true_all_cutoff = y_true_all_cutoff + metrics_summary['y_true_all']
                        y_score_all_cutoff = y_score_all_cutoff + metrics_summary['y_score_all']
            precis[prot_name + '_' + prot_chain] = precis_cutoff
            recall[prot_name + '_' + prot_chain] = recall_cutoff
            f1[prot_name + '_' + prot_chain] = f1_cutoff
            fb[prot_name + '_' + prot_chain] = fb_cutoff
            mcc[prot_name + '_' + prot_chain] = mcc_cutoff
            acc[prot_name + '_' + prot_chain] = acc_cutoff
            h_loss[prot_name + '_' + prot_chain] = h_loss_cutoff
            bal_acc[prot_name + '_' + prot_chain] = bal_acc_cutoff
            jaccard[prot_name + '_' + prot_chain] = jaccard_cutoff
            zo_loss[prot_name + '_' + prot_chain] = zo_loss_cutoff
            tp[prot_name + '_' + prot_chain] = tp_cutoff
            fp[prot_name + '_' + prot_chain] = fp_cutoff
            tn[prot_name + '_' + prot_chain] = tn_cutoff
            fn[prot_name + '_' + prot_chain] = fn_cutoff
            ppv_custom[prot_name + '_' + prot_chain] = ppv_cutoff
            auc[prot_name + '_' + prot_chain] = auc_cutoff
            ap[prot_name + '_' + prot_chain] = ap_cutoff

            roc_fpr[prot_name + '_' + prot_chain] = roc_fpr_cutoff
            roc_tpr[prot_name + '_' + prot_chain] = roc_tpr_cutoff
            roc_thres[prot_name + '_' + prot_chain] = roc_thres_cutoff
            roc_fpr_custom[prot_name + '_' + prot_chain] = roc_fpr_custom_cutoff
            roc_tpr_custom[prot_name + '_' + prot_chain] = roc_tpr_custom_cutoff

            pr_precs[prot_name + '_' + prot_chain] = pr_precis_cutoff
            pr_recall[prot_name + '_' + prot_chain] = pr_recall_cutoff
            pr_thres[prot_name + '_' + prot_chain] = pr_thres_cutoff
            pr_precs_custom[prot_name + '_' + prot_chain] = pr_precis_custom_cutoff
            pr_recall_custom[prot_name + '_' + prot_chain] = pr_recall_custom_cutoff

            y_true_all = y_true_all + y_true_all_cutoff
            y_score_all = y_score_all + y_score_all_cutoff
        # print(y_true_all)
        # print(len(y_true_all))
        if sv_fp:
            with open(sv_fp + tool + '_precis.json', 'w') as filep:
                json.dump(precis, filep)
            with open(sv_fp + tool + '_recall.json', 'w') as filep:
                json.dump(recall, filep)
            with open(sv_fp + tool + '_f1.json', 'w') as filep:
                json.dump(f1, filep)
            with open(sv_fp + tool + '_fb.json', 'w') as filep:
                json.dump(fb, filep)
            with open(sv_fp + tool + '_mcc.json', 'w') as filep:
                json.dump(mcc, filep)
            with open(sv_fp + tool + '_acc.json', 'w') as filep:
                json.dump(acc, filep)
            with open(sv_fp + tool + '_h_loss.json', 'w') as filep:
                json.dump(h_loss, filep)
            with open(sv_fp + tool + '_bal_acc.json', 'w') as filep:
                json.dump(bal_acc, filep)
            with open(sv_fp + tool + '_jaccard.json', 'w') as filep:
                json.dump(jaccard, filep)
            with open(sv_fp + tool + '_zo_loss.json', 'w') as filep:
                json.dump(zo_loss, filep)
            with open(sv_fp + tool + '_tp.json', 'w') as filep:
                print(tp)
                json.dump(tp, filep)
            with open(sv_fp + tool + '_fp.json', 'w') as filep:
                json.dump(fp, filep)
            with open(sv_fp + tool + '_tn.json', 'w') as filep:
                json.dump(tn, filep)
            with open(sv_fp + tool + '_fn.json', 'w') as filep:
                json.dump(fn, filep)
            with open(sv_fp + tool + '_ppv_custom.json', 'w') as filep:
                json.dump(ppv_custom, filep)
            with open(sv_fp + tool + '_auc.json', 'w') as filep:
                json.dump(auc, filep)
            with open(sv_fp + tool + '_ap.json', 'w') as filep:
                json.dump(ap, filep)

            self.pfwriter.generic(pd.DataFrame(y_true_all), sv_fpn=sv_fp + tool + '_y_true_all.txt')
            self.pfwriter.generic(pd.DataFrame(y_score_all), sv_fpn=sv_fp + tool + '_y_score_all.txt')
            # with open(sv_fp + tool + '_y_true_all.json', 'w') as filep:
            #     json.dump(y_true_all, filep)
            # with open(sv_fp + tool + '_y_score_all.json', 'w') as filep:
            #     json.dump(y_score_all, filep)

            with open(sv_fp + tool + '_pr_precs.json', 'w') as filep:
                json.dump(pr_precs, filep)
            with open(sv_fp + tool + '_pr_recall.json', 'w') as filep:
                json.dump(pr_recall, filep)
            with open(sv_fp + tool + '_pr_thres.json', 'w') as filep:
                json.dump(pr_thres, filep)
            with open(sv_fp + tool + '_pr_precs_custom.json', 'w') as filep:
                json.dump(pr_precs_custom, filep)
            with open(sv_fp + tool + '_pr_recall_custom.json', 'w') as filep:
                json.dump(pr_recall_custom, filep)

            with open(sv_fp + tool + '_roc_fpr.json', 'w') as filep:
                json.dump(roc_fpr, filep)
            with open(sv_fp + tool + '_roc_tpr.json', 'w') as filep:
                json.dump(roc_tpr, filep)
            with open(sv_fp + tool + '_roc_thres.json', 'w') as filep:
                json.dump(roc_thres, filep)
            with open(sv_fp + tool + '_roc_fpr_custom.json', 'w') as filep:
                json.dump(roc_fpr_custom, filep)
            with open(sv_fp + tool + '_roc_tpr_custom.json', 'w') as filep:
                json.dump(roc_tpr_custom, filep)
        return 'Finished'

    def inspect(self, df_fpn):
        matrix = self.pfreader.generic(df_fpn=df_fpn)
        # acc_matrix = acc_matrix.reset_index(inplace=False, drop=True)
        # print(matrix.ix[(matrix==0).all(axis=1), :] )
        # matrix = matrix.loc[~(matrix == 0).all(axis=1),:]
        # matrix = matrix.ix[~(matrix == 0).all(axis=1),:]
        return matrix.mean().round(5)

    def pos_sgl_list(
            self,
            prot_name,
            seq_chain,
            pdb_fp,
            topo_fp,
            xml_fp,
            fasta_fp,
            segment,
    ):
        import tmkit as tmk
        pdbtm_seg, pred_seg = tmk.topo.cepdbtm(
            pdb_fp=pdb_fp,
            prot_name=prot_name,
            seq_chain=seq_chain,
            file_chain=chainname().chain(seq_chain),
            topo_fp=topo_fp,
            xml_fp=xml_fp,
            fasta_fp=fasta_fp,
        )
        # print(pdbtm_seg, pred_seg)
        pdbtm_pos_cyto = tmk.seq.pos_seg_list_single(pdbtm_seg['cyto_lower'], pdbtm_seg['cyto_upper'])
        pdbtm_pos_tmh = tmk.seq.pos_seg_list_single(pdbtm_seg['tmh_lower'], pdbtm_seg['tmh_upper'])
        pdbtm_pos_extra = tmk.seq.pos_seg_list_single(pdbtm_seg['extra_lower'], pdbtm_seg['extra_upper'])
        pred_pos_cyto = tmk.seq.pos_seg_list_single(pred_seg['cyto_lower'], pred_seg['cyto_upper'])
        pred_pos_tmh = tmk.seq.pos_seg_list_single(pred_seg['tmh_lower'], pred_seg['tmh_upper'])
        pred_pos_extra = tmk.seq.pos_seg_list_single(pred_seg['extra_lower'], pred_seg['extra_upper'])
        pdbtm_pos = pdbtm_pos_cyto + pdbtm_pos_tmh + pdbtm_pos_extra
        pred_pos = pred_pos_cyto + pred_pos_tmh + pred_pos_extra
        if segment == 'pdbtm_cyto':
            pos_single_list = pdbtm_pos_cyto
            len_seg = fealength().segment(
                seg_lower=pdbtm_seg['cyto_lower'],
                seg_upper=pdbtm_seg['cyto_upper'],
            )
        elif segment == 'pdbtm_tmh':
            pos_single_list = pdbtm_pos_tmh
            len_seg = fealength().segment(
                seg_lower=pdbtm_seg['tmh_lower'],
                seg_upper=pdbtm_seg['tmh_upper'],
            )
        elif segment == 'pdbtm_extra':
            pos_single_list = pdbtm_pos_extra
            len_seg = fealength().segment(
                seg_lower=pdbtm_seg['extra_lower'],
                seg_upper=pdbtm_seg['extra_upper'],
            )
        elif segment == 'phobius_cyto':
            pos_single_list = pred_pos_cyto
            len_seg = fealength().segment(
                seg_lower=pred_seg['cyto_lower'],
                seg_upper=pred_seg['cyto_upper'],
            )
        elif segment == 'phobius_tmh':
            pos_single_list = pred_pos_tmh
            len_seg = fealength().segment(
                seg_lower=pred_seg['tmh_lower'],
                seg_upper=pred_seg['tmh_upper'],
            )
        elif segment == 'phobius_extra':
            pos_single_list = pred_pos_extra
            len_seg = fealength().segment(
                seg_lower=pred_seg['extra_lower'],
                seg_upper=pred_seg['extra_upper'],
            )
        elif segment == 'pdbtm_combined':
            pos_single_list = pdbtm_pos
            part_cyto = fealength().segment(
                seg_lower=pdbtm_seg['cyto_lower'],
                seg_upper=pdbtm_seg['cyto_upper'],
            )
            part_tmh = fealength().segment(
                seg_lower=pdbtm_seg['tmh_lower'],
                seg_upper=pdbtm_seg['tmh_upper'],
            )
            part_extra = fealength().segment(
                seg_lower=pdbtm_seg['extra_lower'],
                seg_upper=pdbtm_seg['extra_upper'],
            )
            len_seg = part_cyto + part_tmh + part_extra
        else:
            pos_single_list = pred_pos
            part_cyto = fealength().segment(
                seg_lower=pred_seg['cyto_lower'],
                seg_upper=pred_seg['cyto_upper'],
            )
            part_tmh = fealength().segment(
                seg_lower=pred_seg['tmh_lower'],
                seg_upper=pred_seg['tmh_upper'],
            )
            part_extra = fealength().segment(
                seg_lower=pred_seg['extra_lower'],
                seg_upper=pred_seg['extra_upper'],
            )
            len_seg = part_cyto + part_tmh + part_extra
        return pos_single_list, len_seg


if __name__ == "__main__":
    from pypropel.path import to

    p = Dispatcher()

    # pos_sgl_list, len_seg = p.pos_sgl_list(
    #     prot_name='1aij',
    #     seq_chain='L',
    #     pdb_fp=to('data/pdb/pdbtm/'),
    #     topo_fp=to('data/phobius/'),
    #     xml_fp=to('data/xml/'),
    #     fasta_fp=to('data/fasta/'),
    #     segment='pdbtm_tmh',
    # )

    p.segment(
        prot_df=pd.DataFrame({
            'prot': ['1aij', ],
            'chain': ['L', ],
        }),
        dist_fp=to('data/pdb/complex/pdbtm/'),
        dist_limit=6.,
        tool_fp=to('data/isite/deeptminter/'),
        tool='tma300',
        pdb_fp=to('data/pdb/pdbtm/'),
        topo_fp=to('data/phobius/'),
        xml_fp=to('data/xml/'),
        fasta_fp=to('data/fasta/'),
        segment='pdbtm_tmh',
        sort=1,
        sv_fp='./',
    )

    # # /*** inspect ***/
    # print(p.inspect(
    #     df_fpn=to('data/al/prediction/ppisite/') + tools[0] + '/' + datasets[0] + '/' + topologies[0] + '/' + regions[0] + '/'
    #            + 'mbpred_fn.txt',
    #     # df_fpn=to('data/al/prediction/ppisite/') + tools[0] + '/' + meet + '/' + datasets[0] + '/' + topologies[0] + '/' + regions[0] + '/' + models[0] + '/tma300_precis.txt'
    # ))