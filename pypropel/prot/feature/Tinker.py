__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import time


class tinker:

    def __init__(self, ):
        pass

    def fill(self, num, list_2d):
        start_time = time.time()
        list_2d_ = list_2d
        len_list_2d = len(list_2d_)
        for i in range(len_list_2d):
            list_2d_[i] = list_2d_[i] + [0 for _ in range(num)]
        end_time = time.time()
        print('------> tinker fill feature: {time}s.'.format(time=end_time - start_time))
        return list_2d_