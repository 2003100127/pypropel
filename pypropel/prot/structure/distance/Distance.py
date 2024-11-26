__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import os
import sys
sys.path.append(os.path.dirname(os.getcwd()) + '/')
from abc import ABCMeta, abstractmethod
# from Bio.PDB.Polypeptide import three_to_one
from pypropel.util.Console import Console
console = Console()

three_to_one = {
    'CYS': 'C',
    'ASP': 'D',
    'SER': 'S',
    'GLN': 'Q',
    'LYS': 'K',
    'ILE': 'I',
    'PRO': 'P',
    'THR': 'T',
    'PHE': 'F',
    'ASN': 'N',
    'GLY': 'G',
    'HIS': 'H',
    'LEU': 'L',
    'ARG': 'R',
    'TRP': 'W',
    'ALA': 'A',
    'VAL':'V',
    'GLU': 'E',
    'TYR': 'Y',
    'MET': 'M',
}


class distance(metaclass=ABCMeta):

    @abstractmethod
    def calculate(self):
        pass

    def one2one_minimal(
            self,
            chain1,
            chain2,
            verbose: bool = False,
    ):
        """
        Notes
        -----
        It outputs minimal distances of each residue in chain 1 to a residue in chain 2.


        Parameters
        ----------
        chain1
            a biopython-typed chain ob 1
        chain2
            a biopython-typed chain ob 2

        Returns
        -------
        2d array

        """
        console.verbose = verbose
        dist_matrix = []
        count_hetamt_1 = 0
        count_hetamt_2 = 0
        for index_1, residue_1 in enumerate(chain1):
            console.print("==================>residue 1 ID: {}".format(index_1))
            if residue_1.get_id()[0] != ' ':
                count_hetamt_1 = count_hetamt_1 + 1
                continue
            else:
                residue_dist = []
                for index_2, residue_2 in enumerate(chain2):
                    # console.print("=====================>residue 2 ID: {}".format(index_1))
                    if residue_2.get_id()[0] != ' ':
                        count_hetamt_2 = count_hetamt_2 + 1
                        continue
                    else:
                        tmp_atom_dist = []
                        for atom_1 in residue_1:
                            if atom_1.get_name() != 'H':
                                for atom_2 in residue_2:
                                    if atom_2.get_name() != 'H':
                                        calc_dist = residue_1[atom_1.get_name()] - residue_2[atom_2.get_name()]
                                        tmp_atom_dist.append(calc_dist)
                        # print('tmp list: %s' % tmp_atom_dist)
                        min_atom_dist = min(tmp_atom_dist)
                        # print(min_atom_dist)
                        residue_dist.append(min_atom_dist)
                # if min(residue_dist) < 6:
                #     print(min(residue_dist))
                min_residue_dist = min(residue_dist)
                dist_matrix.append([
                    index_1 + 1 - count_hetamt_1,
                    three_to_one(residue_1.get_resname()),
                    residue_1.id[1],
                    min_residue_dist
                ])
        return dist_matrix

    def one2one_all(
            self,
            chain1,
            chain2,
            verbose: bool = False,
    ):
        """

        Notes
        -----
        It outputs the minimum distance of each residue in chain 1 to
        each residue in the chain 2.

        Parameters
        ----------
        chain1
            a biopython-typed chain ob 1
        chain2
            a biopython-typed chain ob 2

        Returns
        -------
        2d array
        """
        console.verbose = verbose
        dist_matrix = []
        count_hetamt_1 = 0
        count_hetamt_2 = 0
        for index_1, residue_1 in enumerate(chain1):
            console.print("==================>residue 1 ID: {}".format(index_1))
            if residue_1.get_id()[0] != ' ':
                # print(residue_1.get_id())
                count_hetamt_1 = count_hetamt_1 + 1
                continue
            else:
                # print(residue_1.get_id())
                for index_2, residue_2 in enumerate(chain2):
                    # console.print("=====================>residue 2 ID: {}".format(index_1))
                    tmp_atom_dist = []
                    if residue_2.get_id()[0] != ' ':
                        count_hetamt_2 = count_hetamt_2 + 1
                        continue
                    else:
                        for atom_1 in residue_1:
                            if atom_1.get_name() != 'H':
                                for atom_2 in residue_2:
                                    if atom_2.get_name() != 'H':
                                        calc_dist = residue_1[atom_1.get_name()] - residue_2[atom_2.get_name()]
                                        tmp_atom_dist.append(calc_dist)
                        # print('tmp list: %s' % tmp_atom_dist)
                        min_dist = min(tmp_atom_dist)
                        # print(min_dist)
                        dist_matrix.append([
                            index_1 + 1 - count_hetamt_1,
                            three_to_one(residue_1.get_resname()),
                            residue_1.id[1],
                            index_2 + 1 - count_hetamt_2,
                            'U' if residue_2.get_resname() == 'UNK' else three_to_one(residue_2.get_resname()),
                            residue_2.id[1],
                            min_dist,
                        ])
        return dist_matrix

    def check(
            self,
            chain1,
            chain2,
            thres=6,
            verbose: bool = False,
    ):
        """
        Each residue in the chain 1 has a minimum distance against all of
        residues in the chain 2.
        It stops the calculations of the minimum distance of each residue
        to each residue in the chain 2 when it detects a minimum distance
        of less than thres, 6 by default.

        Parameters
        ----------
        chain1
            a biopython-typed chain ob 1
        chain2
            a biopython-typed chain ob 2
        thres

        Returns
        -------
        2d array

        """
        console.verbose = verbose
        mark = False
        for index_1, residue_1 in enumerate(chain1):
            console.print("==================>residue 1 ID: {}".format(index_1))
            min_check = []
            if residue_1.get_id()[0] != ' ':
                continue
            else:
                for index_2, residue_2 in enumerate(chain2):
                    tmp_atom_dist = []
                    if residue_2.get_id()[0] != ' ':
                        continue
                    else:
                        for atom_1 in residue_1:
                            if atom_1.get_name() != 'H':
                                for atom_2 in residue_2:
                                    if atom_2.get_name() != 'H':
                                        calc_dist = residue_1[atom_1.get_name()] - residue_2[atom_2.get_name()]
                                        tmp_atom_dist.append(calc_dist)
                        # print('tmp list: %s' % tmp_atom_dist)
                        min_dist = min(tmp_atom_dist)
                        min_check.append(min_dist)
                    if min(min_check) < thres:
                        mark = True
                        console.print("==================>residue {} and residue {} in interaction".format(index_1, index_2))
                        break
                if mark:
                    break
        return mark