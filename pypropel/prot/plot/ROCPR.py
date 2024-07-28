__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import seaborn as sns
import matplotlib.pyplot as plt
from pypropel.prot.plot.MetricFormatter import MetricFormatter


class ROCPR:

    def __init__(self, ):
        self.mformatter = MetricFormatter()
        sns.set(font="Helvetica")
        sns.set_style("ticks")

    def draw(
            self,
            X_fpns,
            Y_fpns,
            x_label,
            y_label,
            title,
            ax,
    ):
        methods = [*X_fpns.keys()]
        print(methods)
        prediction_dict = {}
        for method in methods:
            prediction_dict[method] = {}
            x_arr, y_arr = self.mformatter.roc_json_one_prot(X_fpns[method], Y_fpns[method])
            prediction_dict[method]['x'] = x_arr
            prediction_dict[method]['y'] = y_arr
        print(prediction_dict)
        tableau = [
            # 'cornflowerblue',
            # 'steelblue',
            # 'royalblue',
            # 'darkslateblue',
            # 'cadetblue',
            # 'darkcyan',
            # 'teal',
            # 'seagreen',
            # 'lightslategrey',
            # 'darkslategray',
            # # 'steelblue',
            # 'darkseagreen',
            # 'mediumpurple',
            # 'black',

            'lightcoral',
            'darksalmon',
            'indianred',
            'brown',
            'darkred',
            'coral',
            'orangered',
            'chocolate',
            'sienna',
            'saddlebrown',
            'goldenrod',
            'crimson',
            'black',
        ]
        markers = [
            ".",
            "o",
            "v",
            "^",
            "<",
            ">",
            "s",
            "d",
            "D",
            "X_fpns",
            "h",
            "*",
            'x'
        ]
        for i, (method, prediction) in enumerate(prediction_dict.items()):
            if i > 10:
                ax.plot(
                    prediction['x'],
                    prediction['y'],
                    marker=markers[i],
                    markevery=4,
                    label=method,
                    color=tableau[i],
                    linestyle="-",
                    linewidth=2,
                    alpha=0.75,
                )
            else:
                ax.plot(
                    prediction['x'],
                    prediction['y'],
                    marker=markers[i],
                    markevery=10,
                    label=method,
                    color=tableau[i],
                    linestyle="-",
                    linewidth=2,
                    alpha=0.75,
                )
        ax.spines['top'].set_color('none')
        ax.spines['right'].set_color('none')
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title)
        plt.legend(fontsize=10)



if __name__ == "__main__":
    from pypropel.path import to

    X_fpns = {
        'tma300': to('data/eval/tma300/tma300_roc_fpr_custom.json'),
    }
    Y_fpns = {
        'tma300': to('data/eval/tma300/tma300_roc_tpr_custom.json'),
    }
    p = ROCPR()

    fig, ax = plt.subplots(
        nrows=2,
        ncols=2,
        # figsize=(6, 5),
        figsize=(12, 10),
        sharey='all',
        sharex=False,
    )
    print(p.draw(
        X_fpns,
        Y_fpns,
        x_label='fpr',
        y_label='tpr',
        title='',
        ax=ax[0, 0],
    ))
    print(p.draw(
        X_fpns,
        Y_fpns,
        x_label='fpr',
        y_label='tpr',
        title='',
        ax=ax[0, 1],
    ))
    print(p.draw(
        X_fpns,
        Y_fpns,
        x_label='fpr',
        y_label='tpr',
        title='',
        ax=ax[1, 0],
    ))
    plt.show()