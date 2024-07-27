__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import warnings

warnings.filterwarnings('ignore')
import numpy as np
from sklearn import metrics


class Evaluate:

    def __init__(self, ):
        pass

    def accuracy(self, y_true, y_pred):
        return metrics.accuracy_score(y_true=y_true, y_pred=y_pred)

    def accuracytopk(self, y_true, y_score, k=10):
        return metrics.top_k_accuracy_score(y_true, y_score, k=k)

    def precision(self, y_true, y_pred):
        return metrics.precision_score(y_true=y_true, y_pred=y_pred)

    def recall(self, y_true, y_pred):
        return metrics.recall_score(y_true, y_pred)

    def specificity(self, y_true, y_pred):
        return

    def mcc(self, y_true, y_pred):
        return metrics.matthews_corrcoef(y_true, y_pred)

    def f1score(self, y_true, y_pred):
        return metrics.f1_score(y_true, y_pred)

    def fbscore(self, y_true, y_pred, b=0.5):
        return metrics.fbeta_score(y_true, y_pred, beta=b)

    def roc(self, y_true, y_scores):
        fpr, tpr, thres = metrics.roc_curve(
            y_true=y_true,
            y_score=y_scores,
        )
        return fpr, tpr, thres

    def ppv_custom(self, y_true, y_pred):
        tp, fp, tn, fn = self.cfmatrix1d(y_true, y_pred)
        return tp / (tp + fp)

    def roc_custom(self, y_true, y_scores, step=0.01, number=150):
        """
        ..  @description:
            -------------
            fpr: false positive rate
            tpr: true positive rate

        :param y_true:
        :param y_scores:
        :param step:
        :return:
        """
        y_true_ = np.array(y_true)
        y_scores_ = np.array(y_scores)
        thres = np.linspace(min(y_scores_) + 0.005, max(y_scores_), number)
        # thres = np.linspace(0, 1, number)
        # thres = np.arange(0.01, 1, step)
        fpr = []
        tpr = []
        for i in range(thres.shape[0]):
            shape = y_scores_[y_scores_ >= thres[i]].shape[0]
            offset_shape = y_scores_.shape[0] - shape
            y_pred_ = np.concatenate(
                [np.zeros([shape]) + 1, np.zeros([offset_shape])],
                axis=0
            ).astype(np.int64)
            tp, fp, tn, fn = self.cfmatrix1d(y_true_, y_pred_)
            fpr.append(fp / (tn + fp))
            tpr.append(tp / (tp + fn))
        return fpr, tpr

    def roc_custom_L(self, y_true, y_scores, thres):
        """
        ..  @description:
            -------------
            fpr: false positive rate
            tpr: true positive rate

        :param y_true:
        :param y_scores:
        :param step:
        :return:
        """
        y_true_ = np.array(y_true)
        y_scores_ = np.array(y_scores)
        thres = np.array(thres)
        fpr = []
        tpr = []
        for i in range(thres.shape[0]):
            shape = y_scores_[: int(round(thres[i]))].shape[0]
            offset_shape = y_scores_.shape[0] - shape
            y_pred_ = np.concatenate(
                [np.zeros([shape]) + 1, np.zeros([offset_shape])],
                axis=0
            ).astype(np.int64)
            tp, fp, tn, fn = self.cfmatrix1d(y_true_, y_pred_)
            fpr.append(fp / (tn + fp))
            tpr.append(tp / (tp + fn))
        return fpr, tpr

    def prc(self, y_true, y_scores):
        p, r, thres = metrics.precision_recall_curve(
            y_true=y_true,
            probas_pred=y_scores,
        )
        return p, r, thres

    def prc_custom(self, y_true, y_scores, inf=0, sup=0.5, step=0.01, number=100):
        y_true_ = np.array(y_true)
        y_scores_ = np.array(y_scores)
        # thres = np.arange(inf, sup, step)
        thres = np.linspace(np.min(y_scores_) + 0.005, np.max(y_scores_), number)
        p = []
        r = []
        for i in range(thres.shape[0]):
            shape = y_scores_[y_scores_ >= thres[i]].shape[0]
            offset_shape = y_scores_.shape[0] - shape
            y_pred_ = np.concatenate(
                [np.zeros([shape]) + 1, np.zeros([offset_shape])],
                axis=0
            ).astype(np.int64)
            p.append(self.precision(y_true_, y_pred_))
            r.append(self.recall(y_true_, y_pred_))
        return p, r, thres

    def prc_custom_L(self, y_true, y_scores, thres):
        """
        ..  @description:
            zero_l = np.insert(L / np.flipud(np.linspace(1, 50, num1)), 0, 0).tolist()
            l_entire = np.linspace(L + 1, entire, num2).tolist()
            thres = zero_l + l_entire

        :param y_true:
        :param y_scores:
        :param number:
        :return:
        """
        y_true_ = np.array(y_true)
        y_scores_ = np.array(y_scores)
        thres = np.array(thres)
        p = []
        r = []
        for i in range(thres.shape[0]):
            shape = y_scores_[: int(round(thres[i]))].shape[0]
            offset_shape = y_scores_.shape[0] - shape
            y_pred_ = np.concatenate(
                [np.zeros([shape]) + 1, np.zeros([offset_shape])],
                axis=0
            ).astype(np.int64)
            p.append(self.precision(y_true_, y_pred_))
            r.append(self.recall(y_true_, y_pred_))
        return p, r, thres

    def ap(self, y_true, y_scores):
        return metrics.average_precision_score(y_true=y_true, y_score=y_scores)

    def auc(self, y_true, y_scores):
        # auc__ = []
        # try:
        #     auc__.append(metrics.roc_auc_score(y_true=y_true, y_score=y_scores))
        # except ValueError:
        #     print('only one class in label cannot be calculated. ')
        #     pass
        # if len(auc__) == 0:
        #     return 0
        # else:
        #     return auc__[0]
        if len(np.unique(y_true)) == 1:
            return metrics.accuracy_score(y_true, np.rint(y_scores))
        return metrics.roc_auc_score(y_true, y_scores)

    def fowlkes_mallows_score(self, y_true, y_pred):
        return metrics.fowlkes_mallows_score(y_true, y_pred)

    def norm_mi(self, y_true, y_pred):
        return metrics.cluster.normalized_mutual_info_score(y_true, y_pred)

    def cfmatrix1d(self, y_true, y_pred):
        """
        for 2 classes.
        :param y_true:
        :param y_pred:
        :return:
        """
        tn, fp, fn, tp = metrics.confusion_matrix(y_true=y_true, y_pred=y_pred, labels=[0, 1]).ravel()
        return tp, fp, tn, fn

    def hamming_loss(self, y_true, y_pred):
        return metrics.hamming_loss(y_true=y_true, y_pred=y_pred)

    def balanced_accuracy_score(self, y_true, y_pred):
        return metrics.balanced_accuracy_score(y_true=y_true, y_pred=y_pred)

    def jaccard_score(self, y_true, y_pred):
        return metrics.jaccard_score(y_true=y_true, y_pred=y_pred)

    def zero_one_loss(self, y_true, y_pred):
        return metrics.zero_one_loss(y_true=y_true, y_pred=y_pred)

    def ce(self, y_true, y_softmax_score):
        return metrics.log_loss(y_true=y_true, y_pred=y_softmax_score)

    def all_metrics(self, y_true, y_pred, target_names=None):
        """
            {'cls1': {
                'precision':0.5,
                'recall':1.0,
                'f1-score':0.67,
                'support':1
                },
            'cls2': { ... },
            }
        :param y_true:
        :param y_pred:
        :param target_names:
        :return:
        """
        if target_names is None:
            target_names = ['cls1', 'cls2']
        return metrics.classification_report(y_true, y_pred, target_names=target_names)


if __name__ == "__main__":
    p = Evaluate()
    # y_true = [0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1]
    # y_pred = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    # y_true = [0, 1, 1, 1, 1, 0, 0]
    # y_pred = [0, 1, 0, 1, 1, 0, 0]
    # y_true = np.array([1, 1, 1, 1])
    # y_scores = np.array([1, 0, 0, 0])
    y_true = [0, 0, 1, 1]
    y_scores = [0.1, 0.4, 0.35, 0.8]
    # y_true = [-1, +1, +1, +1, +1, -1, -1]
    # y_pred = [-1, +1, -1, +1, +1, -1, -1]
    # y_true = [1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1]
    # y_pred = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    # names = ['cls1', 'cls2']
    # print(p.auc(y_true, y_scores))
    # print(p.cfmatrix1d(y_true, y_pred))
    print(p.roc(y_true, y_scores))
    # print(p.prc(y_true, y_scores))
    # print(p.prc_custom(y_true, y_scores))
    # print(p.all_metrics(y_true, y_pred, names))
    # print(p.precision(y_true, y_pred))
    # print(p.recall(y_true, y_pred))
    # print(p.mcc(y_true, y_pred))
    # print(p.accuracy(y_true, y_pred))
    # print(p.fbscore(y_true, y_pred))
    # print(p.hamming_loss(y_true, y_pred))
    # print(p.jaccard_similarity_score(y_true, y_pred))
    # print(p.zero_one_loss(y_true, y_pred))
    # print(p.balanced_accuracy_score(y_true, y_pred))