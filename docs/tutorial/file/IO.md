## 1. Get file names

We can get the file names in a specific directory using `pp.io.find_from_folder`.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

df = pp.io.find_from_folder(
    file_path=to('data/ex/xml/'),
    suffix='.xml',
    flag=1,
    sv_fpn=None,
    # sv_fpn=to('data/find.txt'),
)
print(df)
```


:material-note-multiple-outline: Output
``` text
29/07/2024 02:31:43 logger: ======>0. Find file (like "Q86V85"): 1a11
29/07/2024 02:31:43 logger: ======>1. Find file (like "Q86V85"): 1a91
29/07/2024 02:31:43 logger: ======>2. Find file (like "Q86V85"): 1afo
...
29/07/2024 02:31:43 logger: ======>8954. Find file (like "Q86V85"): 8wbx
29/07/2024 02:31:43 logger: ======>8955. Find file (like "Q86V85"): 8wd6
         0
0     1a11
1     1a91
2     1afo
3     1aig
4     1aij
...    ...
8951  8u5b
8952  8wam
8953  8wba
8954  8wbx
8955  8wd6

[8956 rows x 1 columns]
```


!!! info

    **`flag`**: which method used to suit file names. Default value: 1

    * 1 - a general function for finding the prefixes of files

    * 2 - separate protein names and chains from file prefixes, like 1atz_A, PDBTM format

    * 3 - separate protein names and chains from file prefixes, like 1atzA

    * 4 - separate protein names and multiple chains from file prefixes, like 1atz_ABCD

    * 5 - separate protein names and multiple chains from file prefixes with regular expression


## 2. Different and repeated files between 1D lists

To identify different and repeated files between two lists `data/pdbtm_alpha_10.02.2023.txt` and `df[0]`, we can do

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

df_lg = pfreader().generic(to('data/pdbtm_alpha_10.02.2023.txt'), df_sep='\t')
df_lg[1], df_lg[2] = zip(*df_lg[0].apply(lambda x: (x.split("_")[0], x.split("_")[1])))
series_prot_names = pd.Series(df_lg[1].unique())

pds_diff, psd_rept = pp.io.list_diff_unipartite(
    pds_lg=series_prot_names,
    pds_sm=df[0],
    sv_diff_fpn=to('data/diff.txt'),
    sv_rept_fpn=to('data/repeat.txt'),
)
print(pds_diff)
print(psd_rept)
```

:material-note-multiple-outline: Output
``` text
# print(pds_diff)
Series([], dtype: object)

# print(psd_rept)
0       3zmj
1       5mur
2       7f92
3       7zdv
4       2fyn
        ... 
8814    8a1v
8815    7wu9
8816    5iou
8817    3rce
8818    7f1s
Length: 8819, dtype: object
```

In addition, we can find newly added proteins between two versions of PDBTM protein lists.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

df_lg = pfreader().generic(to('data/pdbtm_alpha_10.02.2023.txt'), df_sep='\t')
df_sm = pfreader().generic(to('data/pdbtm_alpha_06.30.2023.txt'), df_sep='\t')

pds_diff, psd_rept = pp.io.list_diff_unipartite(
    pds_lg=df_lg[0],
    pds_sm=df_sm[0],
    sv_diff_fpn=to('data/diff.txt'),
    sv_rept_fpn=to('data/repeat.txt'),
)
print(pds_diff)
print(psd_rept)
```

:material-note-multiple-outline: Output
``` text
0       6s7o_A
1       8ha2_A
2       8j5p_P
3       1izl_A
4       8p3w_H
         ...  
2951    8t6u_A
2952    8hju_B
2953    6zjy_H
2954    4qkm_B
2955    6hjr_D
Length: 2956, dtype: object
0        5va3_A
1        8f4f_T
2        6yto_F
3        4knf_D
4        2xq3_A
          ...  
33371    6nt7_B
33372    6jlo_Z
33373    8h9j_8
33374    6idf_E
33375    5jtg_E
Length: 33376, dtype: object
```

## 3. Different and repeated files between 2D lists

We can get the list of Different and repeated files between 2D lists. It drops duplicates by considering two columns in lists.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

df = pd.DataFrame()
df1 = pd.DataFrame()

df[0], df[1] = zip(*pds_diff.apply(lambda x: (x.split("_")[0], x.split("_")[1])))
df1[0], df1[1] = zip(*df_lg[0].apply(lambda x: (x.split("_")[0], x.split("_")[1])))
print(df)
print(df1)

df_differ, df_repeat = pp.io.list_diff_bipartite(
    pds_lg_1=df[0],
    pds_lg_2=df[1],
    pds_sm_1=df1[0],
    pds_sm_2=df1[1],
    sv_diff_fpn=to('data/diff1.txt'),
    sv_rept_fpn=to('data/repeat1.txt'),
)
print(df_differ)
print(df_repeat)
```

:material-note-multiple-outline: Output
``` text
# df_differ
Empty DataFrame
Columns: []
Index: []

# df_repeat
         0  1
0     8t2v  A
1     6b8h  3
2     8p3s  H
3     7thu  K
4     8bx5  C
...    ... ..
2951  7a24  V
2952  6vam  B
2953  8go3  C
2954  8c29  n
2955  6psn  F

[2956 rows x 2 columns]
```

# 4. move, copy, remove, rename, create

PyPropel also provides a few `move`, `copy`, `remove`, `rename`, `create` operations on a list of IDs rendered in pandas Series format.

```text
move_files
copy_files
delete_files
rename_file_suffix
rename_file_prefix
makedir
```