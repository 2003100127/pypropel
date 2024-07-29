__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from pypropel.util.Reader import Reader as pfreader
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


class Variant:

    def __init__(self, ):
        self.pfreader = pfreader()
        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def mutpred2(
            self,
            fpn,
            sheet_name,
            title='',
    ):
        df = self.pfreader.excel(fpn, sheet_name=sheet_name, header=0)
        print(df.columns)
        fig, ax = plt.subplots()

        fruits = df['ID'].values
        counts = df['MutPred2 score'].values
        bar_labels = []
        for i, score in enumerate(df['MutPred2 score'].values):
            if score >= 0.5:
                if 'Disease-causing' in bar_labels:
                    bar_labels.append('_Disease-causing')
                else:
                    bar_labels.append('Disease-causing')
            else:
                if 'Neutral' in bar_labels:
                    bar_labels.append('_Neutral')
                else:
                    bar_labels.append('Neutral')
        bar_colors = df['MutPred2 score'].apply(lambda x: 'tab:pink' if x >= 0.5 else 'teal')
        # bar_colors = ['tab:red', 'tab:blue', 'tab:red', 'tab:orange']

        ax.bar(
            fruits,
            counts,
            width=0.5,
            label=bar_labels,
            color=bar_colors,
            alpha=0.5
        )

        ax.set_ylabel('Pathogenic score', fontsize=16)
        ax.set_title(title, fontsize=14)
        ax.set_xticks(np.arange(len(fruits)), fruits, rotation=10, ha='right')

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.legend()
        plt.subplots_adjust(
            # top=0.92,
            bottom=0.14,
            # left=0.18,
            # right=0.95,
            # hspace=0.20,
            # wspace=0.15,
        )
        plt.show()
        return


if __name__ == "__main__":
    from pypropel.path import to

    p = Variant()
    print(p.mutpred2(
        fpn=to('data/mutpred2.xlsx'),
        sheet_name='SR24_CtoU', # SR24_AtoI SR24_CtoU
        title='SR24_CtoU', # SR24_AtoI SR24_CtoU
    ))