__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import math
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


class Binning:

    def __init__(self, df, key, num_bins, ax):
        self.df = df
        self.key = key
        self.num_bins = num_bins + 1
        self.df['agent'] = self.df[key].apply(
            lambda x: x
            # lambda x: np.log(x)
        )
        self.values = self.df['agent'].values
        self.charac_min = self.df['agent'].min(axis=0)
        self.charac_max = self.df['agent'].max(axis=0)
        self.floor = self.charac_min
        # self.floor = self.charac_min - 1e-10 # if not, there will be a class 0
        self.ceil = self.charac_max
        print(
            'The floor for the charac.: {}\n'.format(self.floor),
            'The ceil for the charac.: {}\n'.format(self.ceil),
        )
        self.bins = np.linspace(self.floor, self.ceil, num=self.num_bins).tolist()
        print('The bins are: {}'.format(self.bins))
        self.bin_ids = self.getid()
        self.ax = ax

    def getid(self, right=True):
        """
        ..  @example:
            ---------
            x = np.array([0, 10.0, 12.4, 15.5, 20.])
            ids = np.digitize(x, bins, right=True)
            for n in range(x.size):
                 print(bins[inds[n] - 1], "<=", x[n], "<", bins[inds[n]])

        :param right:
        :return:
        """
        self.df['cls'] = np.digitize(
            x=self.values,
            bins=self.bins,
            right=right,
        )
        cls = self.df['cls'].values
        # print('The bin classes are: {}'.format(cls.tolist()))
        bin_cls = set(self.df['cls'])
        # print('The unique bin classes are: {}'.format(bin_cls))
        ids = {}
        for i, v in enumerate(bin_cls):
            ids[self.bins[i]] = self.df.loc[self.df['cls'] == v].index.tolist()
        # print('The ids are: {}', ids)
        # for n in range(self.values.size):
        #     print(self.bins[cls[n] - 1], "<=", self.values[n], "<", self.bins[cls[n]])
        return ids

    def draw(self, ):
        # dbs = []
        ips = []
        for _, id_list in self.bin_ids.items():
            # dbs.append(self.df.loc[id_list, 'num_db'].values.mean())
            # ips.append(self.df.loc[id_list, 'num_ip'].values.mean())
            ips.append(self.df.loc[id_list, self.key].values.mean())
        todos = {
            # 'dbs': dbs,
            'ips': ips,
        }
        tableau = [
            plt.cm.Greys,
            plt.cm.Reds
        ]
        labels = [
            # 'NIP in PPI networks',
            'NIP in complexes',
        ]
        # print(todos)

        x = [*self.bin_ids]
        num_ctsts = [len(n) for n in self.bin_ids.values()]
        # print(num_ctsts)
        x_ = []
        for i in range(len(x)):
            x_.append(str(i+1) + ' (ns: ' + str(num_ctsts[i]) + ')')

        # rect1 = axes[0].bar(
        #     np.arange(len(todos['dbs'])),
        #     height=todos['dbs'],
        #     width=0.2,
        #     color=tableau[0](0.1),
        #     # label=labels[0],
        #     alpha=0.6,
        #     edgecolor='black',
        #     linewidth=3,
        # )
        # axes[0].set_ylabel(labels[0], fontsize=11)
        rect2 = self.ax.bar(
            np.arange(len(todos['ips'])),
            height=todos['ips'],
            width=0.2,
            color=tableau[0](0.1),
            # label=labels[1],
            alpha=0.6,
            edgecolor='black',
            linewidth=3,
        )
        self.ax.set_xticks(np.arange(len(self.bin_ids.keys())))
        self.ax.set_xticklabels(x_, fontsize=7.3)
        # self.ax.set_xlabel('Bin (number of interaction sites)', fontsize=12)
        self.ax.set_xlabel('Bin', fontsize=12)
        self.ax.set_ylabel(self.key, fontsize=11)

        self.ax.spines['right'].set_color('white')
        self.ax.spines['top'].set_color('white')
        # axes[0].spines['right'].set_color('white')
        # axes[0].spines['top'].set_color('white')



if __name__ == "__main__":
    from pypropel.path import to
    from pypropel.util.Reader import Reader as pfreader
    from pypropel.util.Writer import Writer as pfwriter

    df = pfreader().generic(to('data/binning/ex.txt'), header=0)
    print(df)
    # pfwriter().generic(df=df[['prot_name', 'chain_name', 'num_db', 'num_ip']], sv_fpn=to('data/binning/ex.txt'), header=True)


    # bin_ids = p.getid()
    # print(bin_ids)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 4), sharey=False, sharex='all')
    Binning(
        df=df,
        key='num_db',
        num_bins=10,
        ax=axes[0],
    ).draw()
    Binning(
        df=df,
        key='num_ip',
        num_bins=10,
        ax=axes[1],
    ).draw()
    plt.show()
