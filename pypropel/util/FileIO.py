__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import os
import re
import shutil
import pandas as pd
from pypropel.util.Reader import Reader as pfreader
from pypropel.util.Writer import Writer as pfwriter
from pypropel.util.Console import Console


class FileIO:

    def __init__(
            self,
            verbose: bool = True,
    ):
        self.pfreader = pfreader()
        self.pfwriter = pfwriter()
        self.console = Console()
        self.console.verbose = verbose

    def find(
            self,
            file_path,
            suffix,
            sv_fpn=None,
            flag=1,
    ):
        """
        Find all files with a suffix specified

        Parameters
        ----------
        file_path
            file path
        suffix
            file suffix
        flag
            by default 1
            1 - a general function for finding the prefixes of files
            2 - separate protein names and chains from file prefixes, like 1atz_A, PDBTM format
            3 - separate protein names and chains from file prefixes, like 1atzA
            4 - separate protein names and multiple chains from file prefixes, like 1atz_ABCD
            5 - separate protein names and multiple chains from file prefixes with regular expression
        sv_fpn
            path to a file with file names found and to be saved

        Returns
        -------

        """
        arr = []
        for i, file in enumerate(os.listdir(file_path)):
            ### file
            # 1atzA.pdb
            file_prefix = os.path.splitext(file)[0]
            file_suffix = os.path.splitext(file)[1]
            ### file_prefix
            # 1atzA
            ### file_suffix
            # .pdb
            if file_suffix == suffix:
                if flag == 1:
                    self.console.print('======>{}. Find file (like "Q86V85"): {}'.format(i, file_prefix[0:]))
                    arr.append([file_prefix[0:]])
                elif flag == 2:
                    self.console.print('======>{}. Find file (like "1atz_A"): {}'.format(i, file_prefix[0:]))
                    prot_name = file_prefix[0:4]
                    chain = file_prefix[5]
                    arr.append([prot_name, chain])
                elif flag == 3:
                    self.console.print('======>{}. Find file (like "1atzA"): {}'.format(i, file_prefix[0:]))
                    prot_name = file_prefix[0:4]
                    chain = file_prefix[4]
                    arr.append([prot_name, chain])
                elif flag == 4:
                    self.console.print('======>{}. Find file (like "1atzABCD"): {}'.format(i, file_prefix[0:]))
                    prot_name = file_prefix[0:4]
                    chain = file_prefix[4:]
                    arr.append([prot_name, chain])
                elif flag == 5:
                    self.console.print('======>{}. Find file using regular expression: {}'.format(i, file_prefix[0:]))
                    prot_name = re.sub(r'[A-Z]', '', file_prefix)
                    chain = re.sub(r'[a-z0-9]', '', file_prefix)
                    arr.append([prot_name, chain])
        if sv_fpn:
            self.pfwriter.generic(df=arr, sv_fpn=sv_fpn)
        return pd.DataFrame(arr)

    def move(
            self,
            pds_mv,
            mv_from_fp,
            mv_to_fp,
            suffix,
    ):
        """

        Parameters
        ----------
        pds_mv
            a pandas Series with files to be moved
        mv_from_fp
            path to files before moved
        mv_to_fp
            path to files to be moved
        suffix
            file suffix

        Returns
        -------
        """
        file_prefix_arr = self.find(file_path=mv_from_fp, suffix=suffix, flag=1)[0].values.tolist()
        for i, file_name in enumerate(pds_mv):
            if file_name in file_prefix_arr:
                self.console.print('======>{}. File {} is in the destination folder'.format(i, file_name))
                shutil.move(
                    mv_from_fp + file_name + suffix,
                    mv_to_fp,
                )
        return "Finished!"

    def copy(
            self,
            pds_cp,
            cp_from_fp,
            cp_to_fp,
            suffix,
        ):
        """

        Parameters
        ----------
        pds_cp
            a pandas Series with files to be copied
        cp_from_fp
            path to files before copied
        cp_to_fp
            path to files to be copied
        suffix
            file suffix

        Returns
        -------

        """
        file_prefix_arr = self.find(file_path=cp_from_fp, suffix=suffix, flag=1)[0].values.tolist()
        for i, file_name in enumerate(pds_cp):
            if file_name in file_prefix_arr:
                self.console.print('======>{}. File {} is in the destination folder'.format(i, file_name))
                shutil.copyfile(
                    cp_from_fp + file_name + suffix,
                    cp_to_fp + file_name + suffix,
                )
        return "Finished!"

    def delete(
            self,
            pds_del,
            del_fp,
            suffix,
    ):
        """

        Parameters
        ----------
        pds_del
            a pandas Series to be deleted
        del_fp
            path to files to be deleted
        suffix
            file suffix

        Returns
        -------

        """
        file_prefix_arr = self.find(file_path=del_fp, suffix=suffix, flag=1)[0].values.tolist()
        for i, file_name in enumerate(pds_del):
            if file_name in file_prefix_arr:
                self.console.print('======>{}. File {} is in the destination folder'.format(i, file_name))
                os.remove(del_fp + file_name + suffix)
        return "Finished!"

    def makedir(
            self,
            dir,
    ):
        """

        Parameters
        ----------
        dir
            a directory to be created

        Returns
        -------

        """
        if not os.path.exists(dir):
            os.makedirs(dir)
        return 'Finished.'

    def rename_suffix(
            self,
            pds_rename,
            rename_fp,
            old_suffix,
            new_suffix,
    ):
        """

        Parameters
        ----------

        pds_rename
            a pandas Series to be renamed
        rename_fp
            path to files to be renamed
        old_suffix
            file suffix to be renamed
        new_suffix
            file suffix after renamed

        Returns
        -------

        """
        file_prefix_arr = self.find(file_path=rename_fp, suffix=old_suffix, flag=1)[0].values.tolist()
        for i, file_name in enumerate(pds_rename):
            if file_name in file_prefix_arr:
                self.console.print('======>{}. File {} is in the destination folder'.format(i, file_name))
                os.rename(
                    os.path.join(rename_fp, file_name + old_suffix),
                    os.path.join(rename_fp, file_name + new_suffix)
                )
        return "Finished!"

    def rename_prefix(
            self,
            pds_old,
            pds_new,
            rename_fp,
            suffix,
    ):
        """

        Parameters
        ----------
        pds_old
            a pandas Series to be renamed
        pds_new
            a pandas Series after renamed
        rename_fp
            path to files to be renamed
        suffix
            file suffix to be renamed
        Returns
        -------

        """
        file_prefix_arr = self.find(file_path=rename_fp, suffix=suffix, flag=1)[0].values.tolist()
        for i, file_name in enumerate(pds_old.values):
            if file_name in file_prefix_arr:
                self.console.print('======>{}. File {} is in the destination folder'.format(i, file_name))
                os.rename(
                    os.path.join(rename_fp, file_name + suffix),
                    os.path.join(rename_fp, pds_new[i] + suffix),
                )
        return "Finished!"