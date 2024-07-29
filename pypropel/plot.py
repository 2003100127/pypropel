__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pypropel.prot.feature.alignment.Conservation import Conservation as conser
from pypropel.prot.feature.isite.Reader import Reader as isitereader
from pypropel.prot.plot.ROCPR import ROCPR
from pypropel.prot.plot.Binning import Binning
from pypropel.prot.plot.Variant import Variant


sns.set(font="Helvetica")
sns.set_style("ticks")

def conservation(
        method : str,
        conser_fpns : dict,
        s : float = 1.0,
        alpha : float = 0.5,
        cmap : str = 'CMRmap_r',
        sv_fpn : str = "./conser.pdf",
):
    if method == 'jsd':
        conser_met = conser().jsd
    else:
        conser_met = conser().jsd
    if len(conser_fpns) == 1:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 3), layout='constrained')
    else:
        fig, ax = plt.subplots(nrows=len(conser_fpns), ncols=1, figsize=(6, 2*int(len(conser_fpns))), layout='constrained')

    ax = [ax] if len(conser_fpns) == 1 else ax

    conser_dict = {}
    for k, conser_fpn in conser_fpns.items():
        conser_dict[k] = conser_met(conser_fpn)

    for i, (k, df) in enumerate(conser_dict.items()):
        x = df['alignment_col'].values + 1
        print(df['alignment_col'].values + 1)
        print(df['score'].values)
        ax[i].scatter(
            x,
            df['score'].values,
            # width,
            s=s,
            edgecolor=None,
            c=df['score'].values,
            cmap=cmap,
            alpha=alpha,
        )
        # ax.bar_label(rects, fmt="{:.1%}", padding=2)
        ax[i].set_ylabel('Conservation', fontsize=14)
        ax[i].set_title(k, fontsize=14)
        ax[i].set_xticks(x[::int(df.shape[0]/10)], x[::int(df.shape[0]/10)])
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
    plt.savefig(sv_fpn, format="pdf", bbox_inches="tight")
    plt.show()
    return


def isite(
        method : str,
        isite_fpns: dict,
        s: float = 1.0,
        alpha: float = 0.5,
        cmap: str = 'CMRmap_r',
        sv_fpn: str = "./isite.pdf",
):
    if method == 'graphppis':
        isite_met = isitereader().graphppis
    else:
        isite_met = isitereader().graphppis
    if len(isite_fpns) == 1:
        fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 3), layout='constrained')
    else:
        fig, ax = plt.subplots(nrows=len(isite_fpns), ncols=1, figsize=(6, 2*int(len(isite_fpns))), layout='constrained')

    ax = [ax] if len(isite_fpns) == 1 else ax

    isite_dict = {}
    for k, isite_fpn in isite_fpns.items():
        isite_dict[k] = isite_met(isite_fpn, )

    for i, (k, df) in enumerate(isite_dict.items()):
        x = df['index'].values + 1
        print(df['index'].values + 1)
        print(df['pred_prob'].values)
        ax[i].scatter(
            x,
            df['pred_prob'].values,
            # width,
            s=s,
            edgecolor=None,
            c=df['pred_prob'].values,
            cmap=cmap,
            alpha=alpha,
        )
        # ax.bar_label(rects, fmt="{:.1%}", padding=2)
        ax[i].set_ylabel('Interaction score', fontsize=14)
        ax[i].set_title(k, fontsize=14)
        ax[i].set_xticks(x[::int(df.shape[0]/10)], x[::int(df.shape[0]/10)])
        ax[i].spines['right'].set_visible(False)
        ax[i].spines['top'].set_visible(False)
    plt.savefig(sv_fpn, format="pdf", bbox_inches="tight")
    plt.show()
    return


def rocpr(
        X_fpns,
        Y_fpns,
        x_label,
        y_label,
        title,
        ax,
):
    return ROCPR().draw(
        X_fpns=X_fpns,
        Y_fpns=Y_fpns,
        x_label=x_label,
        y_label=y_label,
        title=title,
        ax=ax,
    )


def binning(
        df,
        key,
        num_bins,
        ax,
):
    return Binning(
        df=df,
        key=key,
        num_bins=num_bins,
        ax=ax,
    ).draw()


def mutpred2(
        fpn,
        sheet_name,
        title='',
):
    return Variant().mutpred2(
        fpn=fpn,
        sheet_name=sheet_name,
        title=title,
    )


if __name__ == "__main__":
    from pypropel.path import to

    # print(conservation(
    #     method='jsd',
    #     conser_fpns={
    #         'ATAD2_LOC113841329': to('data/conservation/jsd/SR24_AtoI/ATAD2_LOC113841329.jsd'),
    #         'CAMK1G': to('data/conservation/jsd/SR24_AtoI/CAMK1G.jsd'),
    #         'CYP2W1_LOC101804267': to('data/conservation/jsd/SR24_AtoI/CYP2W1_LOC101804267.jsd'),
    #         'KIF27': to('data/conservation/jsd/SR24_AtoI/KIF27.jsd'),
    #         'KIF27_LOC113841629': to('data/conservation/jsd/SR24_AtoI/KIF27_LOC113841629.jsd'),
    #         'LOC119718710': to('data/conservation/jsd/SR24_AtoI/LOC119718710.jsd'),
    #         'RBBP8NL': to('data/conservation/jsd/SR24_AtoI/RBBP8NL.jsd'),
    #     },
    #     # conser_fpns={
    #     #     'CLEC2B_LOC113845378': to('data/conservation/jsd/SR24_CtoU/CLEC2B_LOC113845378.jsd'),
    #     #     'KIF27_LOC113841629': to('data/conservation/jsd/SR24_CtoU/KIF27_LOC113841629.jsd'),
    #     #     'LOC101804340': to('data/conservation/jsd/SR24_CtoU/LOC101804340.jsd'),
    #     #     'ZDHHC20_LOC101792807': to('data/conservation/jsd/SR24_CtoU/ZDHHC20_LOC101792807.jsd'),
    #     # },
    #     cmap='CMRmap_r',
    #     sv_fpn="./A2I_conser.pdf", # A2I_conser C2U_conser
    # ))

    # print(isite(
    #     method='graphppis',
    #     isite_fpns={
    #         'ATAD2_LOC113841329': to('data/isite/graphppis/SR24_AtoI/ATAD2_LOC113841329.txt'),
    #         'CYP2W1_LOC101804267': to('data/isite/graphppis/SR24_AtoI/CYP2W1_LOC101804267.txt'),
    #         'KIF27': to('data/isite/graphppis/SR24_AtoI/KIF27.txt'),
    #         'KIF27_LOC113841629': to('data/isite/graphppis/SR24_AtoI/KIF27_LOC113841629.txt'),
    #         'LOC119718710': to('data/isite/graphppis/SR24_AtoI/LOC119718710.txt'),
    #         'RBBP8NL': to('data/isite/graphppis/SR24_AtoI/RBBP8NL.txt'),
    #         'CAMK1G': to('data/isite/graphppis/SR24_AtoI/RBBP8NL.txt'),
    #     },
    #     # isite_fpns={
    #     #     'CLEC2B_LOC113845378': to('data/isite/graphppis/SR24_CtoU/CLEC2B_LOC113845378.txt'),
    #     #     'KIF27_LOC113841629': to('data/isite/graphppis/SR24_CtoU/KIF27_LOC113841629.txt'),
    #     #     'LOC101804340': to('data/isite/graphppis/SR24_CtoU/LOC101804340.txt'),
    #     #     'ZDHHC20_LOC101792807': to('data/isite/graphppis/SR24_CtoU/ZDHHC20_LOC101792807.txt'),
    #     # },
    #     cmap='coolwarm',
    #     sv_fpn="./A2I_ppi.pdf",  # A2I_ppi C2U_ppi
    # ))

    ### ++++++++++++++++++++++++++++++++++++++++++++
    # X_fpns = {
    #     'tma300': to('data/eval/tma300/tma300_roc_fpr_custom.json'),
    # }
    # Y_fpns = {
    #     'tma300': to('data/eval/tma300/tma300_roc_tpr_custom.json'),
    # }
    # fig, ax = plt.subplots(
    #     nrows=2,
    #     ncols=2,
    #     # figsize=(6, 5),
    #     figsize=(12, 10),
    #     sharey='all',
    #     sharex=False,
    # )
    # print(rocpr(
    #     X_fpns,
    #     Y_fpns,
    #     x_label='fpr',
    #     y_label='tpr',
    #     title='',
    #     ax=ax[0, 0],
    # ))
    # print(rocpr(
    #     X_fpns,
    #     Y_fpns,
    #     x_label='fpr',
    #     y_label='tpr',
    #     title='',
    #     ax=ax[0, 1],
    # ))
    # print(rocpr(
    #     X_fpns,
    #     Y_fpns,
    #     x_label='fpr',
    #     y_label='tpr',
    #     title='',
    #     ax=ax[1, 0],
    # ))
    # plt.show()


    ### ++++++++++++++++++++++++++++++++++++++++++++
    # from pypropel.util.Reader import Reader as pfreader
    # df = pfreader().generic(to('data/binning/ex.txt'), header=0)
    # print(df)
    #
    # fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 4), sharey=False, sharex='all')
    # binning(
    #     df=df,
    #     key='num_db',
    #     num_bins=10,
    #     ax=axes[0],
    # ).draw()
    # binning(
    #     df=df,
    #     key='num_ip',
    #     num_bins=10,
    #     ax=axes[1],
    # ).draw()
    # plt.show()


    ### ++++++++++++++++++++++++++++++++++++++++++++
    mutpred2(
        fpn=to('data/mutpred2.xlsx'),
        sheet_name='SR24_CtoU',  # SR24_AtoI SR24_CtoU
        title='SR24_CtoU',  # SR24_AtoI SR24_CtoU
    )
