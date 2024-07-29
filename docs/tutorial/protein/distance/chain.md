## 1. Distance calculation

We calculate the distance between residues from two given chains. The `pp.dist.one_vs_one` makes a one-against-one comparison between residues.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

dist_mat = pp.dist.one_vs_one(
    pdb_path1=to('data/pdb/complex/pdbtm/'),
    pdb_name1='1aij',
    file_chain1='',
    seq_chain1='L',
    pdb_path2=to('data/pdb/complex/pdbtm/'),
    pdb_name2='1aij',
    file_chain2='',
    seq_chain2='M',
)
df_dist = pd.DataFrame(dist_mat)
df_dist = df_dist.rename(columns={
    0: 'res_fas_id1',
    1: 'res1',
    2: 'res_pdb_id1',
    3: 'res_fas_id2',
    4: 'res2',
    5: 'res_pdb_id2',
    6: 'dist',
})
print(df_dist)
```

:material-note-multiple-outline: Output
``` shell
       res_fas_id1 res1  res_pdb_id1  res_fas_id2 res2  res_pdb_id2       dist
0                1    A            1            1    A            1  33.363850
1                1    A            1            2    E            2  27.019474
2                1    A            1            3    Y            3  29.652020
3                1    A            1            4    Q            4  28.374273
4                1    A            1            5    N            5  29.460558
...            ...  ...          ...          ...  ...          ...        ...
84576          281    G          281       -27423    W          297  41.042576
84577          281    G          281       -27422    G          298  42.089241
84578          281    G          281       -27421    Q          299  43.358837
84579          281    G          281       -27420    N          300  44.846886
84580          281    G          281       -27419    H          301  46.933315

[84581 rows x 7 columns]
```


## 2. Check interaction

Each residue in the chain 1 has a minimum distance against all of residues in the chain 2. It stops the calculations of the minimum distance of each residue to each residue in the chain 2 when it detects a minimum distance of less than `thres`, 6 by default. `pp.dist.check_chain_complex` checks two chains given.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

pp.dist.check_chain_paired(
    pdb_fp1=to('data/pdb/pdbtm/'),
    pdb_fp2=to('data/pdb/pdbtm/'),
    prot_name1='1aij',
    prot_name2='1aij',
    prot_chain1='L',
    prot_chain2='M',
    thres=6.,
    sv_fp=to('data/pdb/pdbtm/'),
)
```

:material-note-multiple-outline: Output
``` shell
D:\Programming\anaconda3\envs\prot\Lib\site-packages\Bio\PDB\PDBParser.py:395: PDBConstructionWarning: Ignoring unrecognized record 'END' at line 2234
  warnings.warn(
28/07/2024 20:02:04 logger: =========>Protein PDB code 1: 1aij
28/07/2024 20:02:04 logger: =========>Protein PDB chain 1: L
28/07/2024 20:02:04 logger: =========>Protein PDB code 2: 1aij
28/07/2024 20:02:04 logger: =========>Protein PDB chain 2: M
28/07/2024 20:02:04 logger: ==================>residue 1 ID: 0
D:\Programming\anaconda3\envs\prot\Lib\site-packages\Bio\PDB\PDBParser.py:395: PDBConstructionWarning: Ignoring unrecognized record 'END' at line 2406
  warnings.warn(
28/07/2024 20:02:04 logger: ==================>residue 0 and residue 252 in interaction
```

The results are saved to the `1aijL_1aijM.pcheck` file. It has the following content.
``` text
1aij	L
1aij	M
```