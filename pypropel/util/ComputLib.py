__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from typing import Any, Dict, List, Union

import itertools


class ComputLib:

    def combo2x2(self, arr):
        """
        repeated 2x2 combination of elements of an array.

        Parameters
        ----------
        array

        Returns
        -------
            amount of combo2x2: L*L.

        """
        combo = []
        ob = itertools.product(arr, repeat=2)
        for i in ob:
            combo.append(list(i))
        return combo

    def reverseInterval(self, inf_arr, sup_arr, max_limit):
        ri_inf_arr = []
        ri_sup_arr = []
        for id, ele in enumerate(inf_arr):
            if id == 0:
                if ele != 1:
                    ri_inf_arr.append(1)
                    ri_sup_arr.append(ele - 1)
            else:
                ri_inf_arr.append(sup_arr[id-1] + 1)
                ri_sup_arr.append(ele - 1)
            if id + 1 == len(inf_arr):
                if sup_arr[id] < max_limit:
                    ri_inf_arr.append(sup_arr[id] + 1)
                    ri_sup_arr.append(max_limit)
        return ri_inf_arr, ri_sup_arr

    def tactic6(
            self,
            arr_2d: List[List[Union[int, float, str]]],
    ) -> Dict[Union[int, float, str], List[Union[float, str]]]:
        """
        Apply tactic 6 to a 2D list and return the result as a dictionary.
        If the second dimension of a 2D list is 2, then it returns a dictionary
        that consists of simply the key -> value by taking the 1st element
        (from the list at the second dimension, e.g., [15, 48]) as key
        and the 2nd element as the value. If the list at the second dimension
        has a length of more than 2, then the value of the returned dictionary
        will be this list but eliminating its 1st element, because
        the 1st element will be taken as the key. See Examples below.

        Parameters
        ----------
        arr_2d : List[List[Union[int, float, str]]]
            The input 2D array.

        Examples
        --------
        1st example
        >>> arr_2d = [[15, 48], [30, 53], [2, 3]]
        >>> arr_2d
        # output
        >>> {15: 53, 2: 3}

        2nd example
        >>> arr_2d = [[15, 48, 78], [30, 53, 99], [30, 3, 11]]
        >>> arr_2d
        # output
        >>> {15: [48, 78], 30: [3, 11]}

        Returns
        -------
        Dict[Union[int, float, str], List[Union[float, str]]]
            The result dictionary where the first element of each sublist in `arr_2d`
            is the key and the remaining elements are the values.
        """
        result = {}
        len_arr = len(arr_2d[0])
        if len_arr == 2:
            for item in arr_2d:
                result[item[0]] = item[1]
        else:
            for item in arr_2d:
                result[item[0]] = item[1:]
        return result
