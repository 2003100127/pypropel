__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import matplotlib.pyplot as plt
from pypropel.prot.feature.BasicUnit import BasicUnit

import numpy as np

class jsdconservation:

    def __init__(self, ):
        self.bu = BasicUnit()

    def boxplot(self, masp_fpns):
        helix_start_masp1_jsd = self.bu.helix(self.bu.masp1)[0]
        helix_end_masp1_jsd = self.bu.helix(self.bu.masp1)[1]
        sheet_start_masp1_jsd = self.bu.sheet(self.bu.masp1)[0]
        sheet_end_masp1_jsd = self.bu.sheet(self.bu.masp1)[1]
        turn_start_masp1_jsd = self.bu.turn(self.bu.masp1)[0]
        turn_end_masp1_jsd = self.bu.turn(self.bu.masp1)[1]
        masp1_jsd = self.bu.localize(
            helix_start=helix_start_masp1_jsd,
            helix_end=helix_end_masp1_jsd,
            sheet_start=sheet_start_masp1_jsd,
            sheet_end=sheet_end_masp1_jsd,
            turn_start=turn_start_masp1_jsd,
            turn_end=turn_end_masp1_jsd,
            conser_fpn=masp_fpns['masp1_jsd']
        )
        masp1_jsd_helix = masp1_jsd.groupby('ss_cls').get_group('H')[1]
        masp1_jsd_sheet = masp1_jsd.groupby('ss_cls').get_group('B')[1]
        masp1_jsd_other = masp1_jsd.groupby('ss_cls').get_group('Other')[1]
        masp1_jsd_helix_ = masp1_jsd_helix.values.tolist()
        masp1_jsd_sheet_ = masp1_jsd_sheet.values.tolist()
        masp1_jsd_other_ = masp1_jsd_other.values.tolist()
        masp1_jsd_helix__ = masp1_jsd_helix[masp1_jsd_helix.apply(lambda x: x > 0.3)].values.tolist()
        masp1_jsd_sheet__ = masp1_jsd_sheet[masp1_jsd_sheet.apply(lambda x: x > 0.3)].values.tolist()
        masp1_jsd_other__ = masp1_jsd_other[masp1_jsd_other.apply(lambda x: x > 0.3)].values.tolist()
        masp1_jsd_ = [
            masp1_jsd_helix_,
            masp1_jsd_sheet_,
            masp1_jsd_other_
        ]
        masp1_jsd__ = [
            masp1_jsd_helix__,
            masp1_jsd_sheet__,
            masp1_jsd_other__
        ]

        masp1_jsd_len_ = [
            len(masp1_jsd_helix_),
            len(masp1_jsd_sheet_),
            len(masp1_jsd_other_),
        ]
        sum_stat_masp1_ = sum(masp1_jsd_len_)
        masp1_jsd_len__ = [
            len(masp1_jsd_helix__),
            len(masp1_jsd_sheet__),
            len(masp1_jsd_other__),
        ]
        sum_stat_masp1__ = sum(masp1_jsd_len__)

        helix_start_masp2_jsd = self.bu.helix(self.bu.masp2)[0]
        helix_end_masp2_jsd = self.bu.helix(self.bu.masp2)[1]
        sheet_start_masp2_jsd = self.bu.sheet(self.bu.masp2)[0]
        sheet_end_masp2_jsd = self.bu.sheet(self.bu.masp2)[1]
        turn_start_masp2_jsd = self.bu.turn(self.bu.masp2)[0]
        turn_end_masp2_jsd = self.bu.turn(self.bu.masp2)[1]
        masp2_jsd = self.bu.localize(
            helix_start=helix_start_masp2_jsd,
            helix_end=helix_end_masp2_jsd,
            sheet_start=sheet_start_masp2_jsd,
            sheet_end=sheet_end_masp2_jsd,
            turn_start=turn_start_masp2_jsd,
            turn_end=turn_end_masp2_jsd,
            conser_fpn=masp_fpns['masp2_jsd']
        )
        masp2_jsd_helix = masp2_jsd.groupby('ss_cls').get_group(('H'))[1]
        masp2_jsd_sheet = masp2_jsd.groupby('ss_cls').get_group('B')[1]
        masp2_jsd_turn = masp2_jsd.groupby('ss_cls').get_group('T')[1]
        masp2_jsd_other = masp2_jsd.groupby('ss_cls').get_group('Other')[1]
        masp2_jsd_helix_ = masp2_jsd_helix.values.tolist()
        masp2_jsd_sheet_ = masp2_jsd_sheet.values.tolist()
        masp2_jsd_turn_ = masp2_jsd_turn.values.tolist()
        masp2_jsd_other_ = masp2_jsd_other.values.tolist()
        masp2_jsd_helix__ = masp2_jsd_helix[masp2_jsd_helix.apply(lambda x: x > 0.3)].values.tolist()
        masp2_jsd_sheet__ = masp2_jsd_sheet[masp2_jsd_sheet.apply(lambda x: x > 0.3)].values.tolist()
        masp2_jsd_turn__ = masp2_jsd_turn[masp2_jsd_turn.apply(lambda x: x > 0.3)].values.tolist()
        masp2_jsd_other__ = masp2_jsd_other[masp2_jsd_other.apply(lambda x: x > 0.3)].values.tolist()
        masp2_jsd_ = [
            masp2_jsd_helix_,
            masp2_jsd_sheet_,
            masp2_jsd_turn_,
            masp2_jsd_other_
        ]
        masp2_jsd__ = [
            masp2_jsd_helix__,
            masp2_jsd_sheet__,
            masp2_jsd_turn__,
            masp2_jsd_other__
        ]

        masp2_jsd_len_ = [
            len(masp2_jsd_helix_),
            len(masp2_jsd_sheet_),
            len(masp2_jsd_turn_),
            len(masp2_jsd_other_)
        ]
        sum_stat_masp2_ = sum(masp2_jsd_len_)
        masp2_jsd_len__ = [
            len(masp2_jsd_helix__),
            len(masp2_jsd_sheet__),
            len(masp2_jsd_turn__),
            len(masp2_jsd_other__),
        ]
        sum_stat_masp2__ = sum(masp2_jsd_len__)

        x_tick_lables_masp1 = [
            r'$\alpha$-helix',
            r'$\beta$-sheet',
            'unordered'
        ]
        x_tick_lables_masp2 = [
            r'$\alpha$-helix',
            r'$\beta$-sheet',
            r'$\beta$-turn',
            'unordered'
        ]
        filler_palette_masp1 = [
            'skyblue',
            'gray',
            'deeppink',
        ]
        filler_palette_masp2 = [
            'skyblue',
            'gray',
            'sienna',
            'deeppink',
        ]
        x_ticks_masp1 = np.linspace(0, 2, 3)
        x_ticks_masp2 = np.linspace(0, 3, 4)
        meanpointprops = {
            'marker': 'o',
            'markeredgecolor': 'firebrick',
            'markerfacecolor': 'firebrick',
            'markersize': 6
        }
        fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(9, 7), sharex='col', sharey='row')
        fig.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.20, wspace=0.15)
        # ###/*** plot 1-1 ***/
        bplot11 = axes[0, 0].boxplot(
            masp1_jsd_,
            positions=[0, 1, 2],
            showmeans=True,
            meanprops=meanpointprops,
            patch_artist=True
        )
        for i in [0, 1, 2]:
            # x = x_masp1_jsd[i]
            y = masp1_jsd_[i]
            x = np.random.normal(i, 0.05, len(y))
            axes[0, 0].scatter(
                x,
                y,
                color=filler_palette_masp1[i],
                label=str(masp1_jsd_len_[i]) + '\n' + str(round(masp1_jsd_len_[i] / sum_stat_masp1_, 2)) + '%',
                alpha=0.2,
                s=24
            )
        axes[0, 0].yaxis.grid(True)
        axes[0, 0].spines['right'].set_color('none')
        axes[0, 0].spines['top'].set_color('none')
        axes[0, 0].set_ylabel('Conservation')
        axes[0, 0].set_title(r'all values for $\mathtt{[MaSp1]_{12}}$')
        box = axes[0, 0].get_position()
        axes[0, 0].set_position([box.x0, box.y0, box.width, box.height * 0.95])
        handles, labels = axes[0, 0].get_legend_handles_labels()
        axes[0, 0].legend(
            handles,
            labels,
            loc='upper center',
            bbox_to_anchor=(0.5, 1.3),
            ncol=4,
            fontsize=8.5
        )

        # # ###/*** plot 1-2 ***/
        bplot12 = axes[0, 1].boxplot(
            masp2_jsd_,
            positions=[0, 1, 2, 3],
            showmeans=True,
            meanprops=meanpointprops,
            patch_artist=True
        )
        for i in [0, 1, 2, 3]:
            # x = x_masp1_jsd[i]
            y = masp2_jsd_[i]
            x = np.random.normal(i, 0.05, len(y))
            axes[0, 1].scatter(
                x,
                y,
                color=filler_palette_masp2[i],
                label=str(masp2_jsd_len_[i]) + '\n' + str(round(masp2_jsd_len_[i]/sum_stat_masp2_, 2))+'%',
                alpha=0.2,
                s=24
            )
        axes[0, 1].yaxis.grid(True)
        axes[0, 1].spines['right'].set_color('none')
        axes[0, 1].spines['top'].set_color('none')
        axes[0, 1].set_title(r'all values for $\mathtt{[MaSp2]_{32}}$')
        box = axes[0, 1].get_position()
        axes[0, 1].set_position([box.x0, box.y0, box.width, box.height * 0.95])
        handles, labels = axes[0, 1].get_legend_handles_labels()
        axes[0, 1].legend(
            handles,
            labels,
            loc='upper center',
            bbox_to_anchor=(0.5, 1.3),
            ncol=4,
            fontsize=8.5
        )

        # # ###/*** plot 2-1 ***/
        bplot21 = axes[1, 0].boxplot(
            masp1_jsd__,
            positions=[0, 1, 2],
            showmeans=True,
            meanprops=meanpointprops,
            patch_artist=True
        )
        for i in [0, 1, 2,]:
            # x = x_masp1_jsd[i]
            y = masp1_jsd__[i]
            x = np.random.normal(i, 0.1, len(y))
            axes[1, 0].scatter(
                x,
                y,
                color=filler_palette_masp1[i],
                label=str(masp1_jsd_len__[i]) + '\n' + str(round(masp1_jsd_len__[i]/sum_stat_masp1__, 2))+'%',
                alpha=0.2,
                s=24
            )
        axes[1, 0].yaxis.grid(True)
        axes[1, 0].spines['right'].set_color('none')
        axes[1, 0].spines['top'].set_color('none')
        axes[1, 0].set_ylabel('Conservation')
        axes[1, 0].set_xticks(x_ticks_masp1)
        axes[1, 0].set_xticklabels(x_tick_lables_masp1, fontsize=9)
        axes[1, 0].set_title(r'>0.3 values for $\mathtt{[MaSp1]_{12}}$')
        box = axes[1, 0].get_position()
        axes[1, 0].set_position([box.x0, box.y0, box.width, box.height * 0.95])
        handles, labels = axes[1, 0].get_legend_handles_labels()
        axes[1, 0].legend(
            handles,
            labels,
            loc='upper center',
            bbox_to_anchor=(0.5, 1.3),
            ncol=4,
            fontsize=8.5
        )

        # # ###/*** plot 2-2 ***/
        # axes[1, 1].set_xlim(0, 12)
        # axes[1, 1].set_ylim(0.29, 0.45)
        bplot22 = axes[1, 1].boxplot(
            masp2_jsd__,
            positions=[0, 1, 2, 3],
            showmeans=True,

            meanprops=meanpointprops,
            patch_artist=True
        )
        for i in [0, 1, 2, 3]:
            # x = x_masp1_jsd[i]
            y = masp2_jsd__[i]
            x = np.random.normal(i, 0.05, len(y))
            axes[1, 1].scatter(
                x,
                y,
                color=filler_palette_masp2[i],
                label=str(masp2_jsd_len__[i]) + '\n' + str(round(masp2_jsd_len__[i]/sum_stat_masp2__, 2))+'%',
                alpha=0.2,
                s=24
            )
        axes[1, 1].yaxis.grid(True)
        axes[1, 1].spines['right'].set_color('none')
        axes[1, 1].spines['top'].set_color('none')
        axes[1, 1].set_title(r'>0.3 values for $\mathtt{[MaSp2]_{32}}$')
        box = axes[1, 1].get_position()
        axes[1, 1].set_position([box.x0, box.y0, box.width, box.height * 0.95])
        handles, labels = axes[1, 1].get_legend_handles_labels()
        axes[1, 1].legend(
            handles,
            labels,
            loc='upper center',
            bbox_to_anchor=(0.5, 1.3),
            ncol=4,
            fontsize=8.5
        )

        axes[1, 1].set_xticks(x_ticks_masp2)
        axes[1, 1].set_xticklabels(x_tick_lables_masp2, fontsize=9)

        for bplot in (bplot11, bplot21):
            for patch, color in zip(bplot['boxes'], filler_palette_masp1):
                patch.set(color='black', linewidth=1.1)
                patch.set(facecolor=color, alpha=0.3)
            for whisker in bplot['whiskers']:
                whisker.set(color='black', linewidth=1)
            for cap in bplot['caps']:
                cap.set(color='black', linewidth=2)
            for median in bplot['medians']:
                median.set(color='black', linewidth=2)
            for flier in bplot['fliers']:
                flier.set(marker='o', color='y', alpha=0.5)

        for bplot in (bplot12, bplot22):
            for patch, color in zip(bplot['boxes'], filler_palette_masp2):
                patch.set(color='black', linewidth=1.1)
                patch.set(facecolor=color, alpha=0.3)
            for whisker in bplot['whiskers']:
                whisker.set(color='black', linewidth=1)
            for cap in bplot['caps']:
                cap.set(color='black', linewidth=2)
            for median in bplot['medians']:
                median.set(color='black', linewidth=2)
            for flier in bplot['fliers']:
                flier.set(marker='o', color='y', alpha=0.5)
        fig.tight_layout()
        plt.show()
        return 0


if __name__ == "__main__":
    from pypropel.path import to

    p = jsdconservation()
    masp_fpns = {
        'masp1_jsd': to('data/jsd/MaSp1.jsd'),
        'masp2_jsd': to('data/jsd/MaSp2.jsd'),
    }
    print(p.boxplot(masp_fpns))
