__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

from typing import Tuple

import pandas as pd
from pypropel.util.FileIO import FileIO
from pypropel.util.ListDiffer import ListDiffer


def list_diff_unipartite(
        pds_lg : pd.Series,
        pds_sm : pd.Series,
        sv_diff_fpn : str = None,
        sv_rept_fpn : str = None,
) -> Tuple[pd.Series, pd.Series]:
    """

    Parameters
    ----------
    pds_lg
        the 1st longer pandas Series
    pds_sm
        the 2nd shorter pandas Series
    sv_diff_fpn
        file contains different content between two lists
    sv_rept_fpn
        file contains repeated content between two lists

    Returns
    -------

    """
    pds_differ, pds_repeat = ListDiffer().unipartite(
        pds_lg=pds_lg,
        pds_sm=pds_sm,
        sv_diff_fpn=sv_diff_fpn,
        sv_rept_fpn=sv_rept_fpn,
    )
    return pds_differ, pds_repeat


def list_diff_bipartite(
        pds_lg_1 : pd.Series,
        pds_lg_2 : pd.Series,
        pds_sm_1 : pd.Series,
        pds_sm_2 : pd.Series,
        sv_diff_fpn : str = None,
        sv_rept_fpn : str = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """

    Parameters
    ----------
    pds_lg_1
        the 1st longer pandas Series contains the 1st column
    pds_lg_2
        the 2nd shorter pandas Series contains the 2nd column
    pds_sm_1
        the 1st longer pandas Series contains the 1st column
    pds_sm_2
        the 2nd shorter pandas Series contains the 2nd column
    sv_diff_fpn
        file contains different content between two lists
    sv_rept_fpn
        file contains repeated content between two lists

    Returns
    -------

    """
    df_differ, df_repeat = ListDiffer().bipartite(
        pds_lg_1=pds_lg_1,
        pds_lg_2=pds_lg_2,
        pds_sm_1=pds_sm_1,
        pds_sm_2=pds_sm_2,
        sv_diff_fpn=sv_diff_fpn,
        sv_rept_fpn=sv_rept_fpn,
    )
    return df_differ, df_repeat


def find_from_folder(
        file_path : str,
        suffix : str,
        flag : int = 1,
        sv_fpn : str = None,
) -> pd.DataFrame:
    """

    Parameters
    ----------
    file_path
        file path
    suffix
        file suffix
    flag
        which method used to suit file names
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
    df = FileIO().find(
        file_path=file_path,
        suffix=suffix,
        flag=flag,
        sv_fpn=sv_fpn,
    )
    return df


def move_files(
        pds_mv : pd.Series,
        mv_from_fp : str,
        mv_to_fp : str,
        suffix : str,
) -> str:
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
    FileIO().move(
        pds_mv=pds_mv,
        mv_from_fp=mv_from_fp,
        mv_to_fp=mv_to_fp,
        suffix=suffix,
    )
    return "Finished!"


def copy_files(
        pds_cp : pd.Series,
        cp_from_fp : str,
        cp_to_fp : str,
        suffix : str,
) -> str:
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
    FileIO().copy(
        pds_cp=pds_cp,
        cp_from_fp=cp_from_fp,
        cp_to_fp=cp_to_fp,
        suffix=suffix,
    )
    return "Finished!"


def delete_files(
        pds_del : pd.Series,
        del_fp : str,
        suffix : str,
) -> str:
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
    FileIO().delete(
        pds_del=pds_del,
        del_fp=del_fp,
        suffix=suffix,
    )
    return "Finished!"


def rename_file_suffix(
        pds_rename : pd.Series,
        rename_fp : str,
        old_suffix : str,
        new_suffix : str,
) -> str:
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
    FileIO().rename_suffix(
        pds_rename=pds_rename,
        rename_fp=rename_fp,
        old_suffix=old_suffix,
        new_suffix=new_suffix,
    )
    return "Finished!"


def makedir(
        dir : str,
) -> str:
    """

    Parameters
    ----------
    dir
        a directory to be created

    Returns
    -------

    """
    FileIO().makedir(dir=dir)
    return "Finished!"


def rename_file_prefix(
        pds_old : pd.Series,
        pds_new : pd.Series,
        rename_fp : str,
        suffix : str,
) -> str:
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
    FileIO().rename_prefix(
        pds_old=pds_old,
        pds_new=pds_new,
        rename_fp=rename_fp,
        suffix=suffix,
    )
    return "Finished!"


if __name__ == "__main__":
    from pypropel.util.Reader import Reader as pfreader
    from pypropel.path import to

    # df = find_from_folder(
    #     file_path=to('data/ex/xml/'),
    #     suffix='.xml',
    #     flag=1,
    #     sv_fpn=None,
    #     # sv_fpn=to('data/find.txt'),
    # )
    # print(df)

    # df_lg = pfreader().generic(to('data/pdbtm_alpha_10.02.2023.txt'), df_sep='\t')
    # df_lg[1], df_lg[2] = zip(*df_lg[0].apply(lambda x: (x.split("_")[0], x.split("_")[1])))
    # print(pd.Series(df_lg[1].unique()))
    # series_prot_names = pd.Series(df_lg[1].unique())
    #
    # pds_diff, psd_rept = list_diff_unipartite(
    #     pds_lg=series_prot_names,
    #     pds_sm=df[0],
    #     sv_diff_fpn=to('data/diff.txt'),
    #     sv_rept_fpn=to('data/repeat.txt'),
    # )
    # print(pds_diff)
    # print(psd_rept)

    # print(move_files(
    #     pds_mv=df_lg[1].unique(),
    #     mv_from_fp=to('data/ex/pdbtm/'),
    #     mv_to_fp=to('data/ex/tmp/'),
    #     suffix='.pdb',
    # ))

    # print(copy_files(
    #     pds_cp=df_lg[1].unique(),
    #     cp_from_fp=to('data/ex/pdbtm/'),
    #     cp_to_fp=to('data/ex/tmp/'),
    #     suffix='.pdb',
    # ))

    # print(p.delete(
    #     pds_del=df_lg[1].unique(),
    #     del_fp=to('data/ex/tmp/'),
    #     suffix='.pdb',
    # ))

    # print(p.rename_suffix(
    #     pds_rename=df_lg[1].unique(),
    #     rename_fp=to('data/ex/tmp/'),
    #     old_suffix='.pdb',
    #     new_suffix='.ent',
    # ))

    # print(p.rename_prefix(
    #     pds_old=df_lg[1],
    #     pds_new=df_lg[1].apply(lambda x: x + 'hello'),
    #     rename_fp=to('data/ex/tmp/'),
    #     suffix='.pdb',
    # ))


    # # #/*** unipartite ***/
    # df_lg = pfreader().generic(to('data/pdbtm_alpha_10.02.2023.txt'), df_sep='\t')
    # df_sm = pfreader().generic(to('data/pdbtm_alpha_06.30.2023.txt'), df_sep='\t')
    #
    # pds_diff, psd_rept = list_diff_unipartite(
    #     pds_lg=df_lg[0],
    #     pds_sm=df_sm[0],
    #     sv_diff_fpn=to('data/diff.txt'),
    #     sv_rept_fpn=to('data/repeat.txt'),
    # )
    # # print(pds_diff)
    # # print(psd_rept)


    # # #/*** bipartite ***/
    # df = pd.DataFrame()
    # df1 = pd.DataFrame()
    #
    # df[0], df[1] = zip(*pds_diff.apply(lambda x: (x.split("_")[0], x.split("_")[1])))
    # df1[0], df1[1] = zip(*df_lg[0].apply(lambda x: (x.split("_")[0], x.split("_")[1])))
    # print(df)
    # print(df1)
    #
    # df_differ, df_repeat = list_diff_bipartite(
    #     pds_lg_1=df[0],
    #     pds_lg_2=df[1],
    #     pds_sm_1=df1[0],
    #     pds_sm_2=df1[1],
    #     sv_diff_fpn=to('data/diff1.txt'),
    #     sv_rept_fpn=to('data/repeat1.txt'),
    # )
    # print(df_differ)
    # print(df_repeat)