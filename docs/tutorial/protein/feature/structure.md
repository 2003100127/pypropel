Features extracted based on 3D protein structures (or PDB files) are informative to characterise proteins. After loads of protein structures predicted by AlphaFold, AlphaFold2, and AlphaFold3 are made available online, the functions embedded in computational tools to access these structural features are necessary.

!!! warning note
    
    PyPropel is endeavouring to make the extraction of more structure-based features available. Currently, we include relative solvent accessibility (RSA) by [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) and encoded sequences by the [3di](https://github.com/steineggerlab/foldseek) technique.


## Relative solvent accessibility

Here, we take six proteins as an example to generate their RSA files using DSSP. It can be done as follows.

:material-language-python: Python
``` py linenums="1"
import pandas as pd

prot_df = pd.DataFrame({
    'prot': ['3pux', '3rko', '3udc', '3vr8', '4kjs', '4pi2', ],
    'chain': ['G', 'A', 'A', 'D', 'A', 'C', ],
})
```

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

for i in prot_df.index:
    print('No.{}: protein: {} chain: {}'.format(i + 1, prot_df.loc[i, 'prot'], prot_df.loc[i, 'chain']))
    dssp_rsa_run(
        prot_name=prot_df.loc[i, 'prot'],
        prot_chain=prot_df.loc[i, 'chain'],
        pdb_fp='data/pdb/pdbtm/',
        sv_fp='data/rsa/',
    )
```

Then, we can access the RSA file of each protein.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp
import pandas as pd

prot_df = pd.DataFrame({
    'prot': ['3pux', '3rko', '3udc', '3vr8', '4kjs', '4pi2', ],
    'chain': ['G', 'A', 'A', 'D', 'A', 'C', ],
})
for i in prot_df.index:
    print('No.{}: protein: {} chain: {}'.format(i + 1, prot_df.loc[i, 'prot'], prot_df.loc[i, 'chain']))
    df_rsa = dssp_rsa_access(
        prot_name=prot_df.loc[i, 'prot'],
        prot_chain=prot_df.loc[i, 'chain'],
        rsa_fp=to('data/rsa/')
    )
    print(df_rsa)
```

:material-note-multiple-outline: Output
``` text
No.1: protein: 3pux chain: G
     rsa_fas_id  rsa_pdb_ids  rsa_prob
0             1            2  1.000000
1             2            3  0.829787
2             3            4  1.000000
3             4            5  0.696970
4             5            6  0.705882
..          ...          ...       ...
288         289          292  1.000000
289         290          293  0.547619
290         291          294  0.394366
291         292          295  1.000000
292         293          296  1.000000

[293 rows x 3 columns]
No.2: protein: 3rko chain: A
    rsa_fas_id  rsa_pdb_ids  rsa_prob
0            1           15  1.000000
1            2           16  0.695431
2            3           17  0.528302
3            4           18  0.668639
4            5           19  0.690355
..         ...          ...       ...
90          91          122  0.018868
91          92          123  0.682927
92          93          124  0.687117
93          94          125  0.731278
94          95          126  0.823944

[95 rows x 3 columns]
No.3: protein: 3udc chain: A
     rsa_fas_id  rsa_pdb_ids  rsa_prob
0             1           13  1.000000
1             2           14  0.760736
2             3           15  0.804734
3             4           16  0.556098
4             5           17  0.650943
..          ...          ...       ...
262         263          275  0.974522
263         264          276  0.969543
264         265          277  0.917073
265         266          278  0.911290
266         267          279  1.000000

[267 rows x 3 columns]
No.4: protein: 3vr8 chain: D
     rsa_fas_id  rsa_pdb_ids  rsa_prob
0             1           28  1.000000
1             2           29  0.861538
2             3           30  0.915094
3             4           31  0.933962
4             5           32  0.985915
..          ...          ...       ...
124         125          152  0.468085
125         126          153  0.570423
126         127          154  0.863436
127         128          155  0.804124
128         129          156  1.000000

[129 rows x 3 columns]
No.5: protein: 4kjs chain: A
     rsa_fas_id  rsa_pdb_ids  rsa_prob
0             1            4  0.840237
1             2            5  0.456853
2             3            6  0.720812
3             4            7  0.698225
4             5            8  0.408537
..          ...          ...       ...
315         316          347  0.380952
316         317          348  0.035533
317         318          349  0.461929
318         319          350  0.798780
319         320          351  0.756098

[320 rows x 3 columns]
No.6: protein: 4pi2 chain: C
     rsa_fas_id  rsa_pdb_ids  rsa_prob
0             1           16  0.938144
1             2           17  0.561538
2             3           18  0.211268
3             4           19  0.676056
4             5           20  0.441718
..          ...          ...       ...
223         224          252  0.556098
224         225          253  0.695122
225         226          254  0.914634
226         227          255  0.732394
227         228          256  0.958763

[228 rows x 3 columns]
```

## 3di encoded sequence

3Di-encoded sequences refer to sequences encoded with 3D structural information of biomolecules, often proteins, to integrate their three-dimensional spatial features into a simplified format suitable for computational analysis, such as machine learning.

This encoding bridges the gap between a protein's primary sequence (linear amino acid sequence) and its three-dimensional conformation (folded structure). The goal of 3Di encoding is to capture spatial relationships, structural motifs, and functional regions that are critical for biological activity.

Includes 3D spatial features derived from the protein's folded structure, typically obtained from crystallography, NMR, or computational modeling (e.g. Angles and torsions: Phi (ϕ), Psi (ψ), and Omega (ω) angles).

Here, we take 3 proteins as an example to generate their 3di.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp
import pandas as pd

prot_df = pd.DataFrame({
    'prot': ['1aij', '1aig', '1xqf', ],
    'chain': ['L', 'L', 'A', ],
})
for i in prot_df.index:
    print('No.{}: protein: {} chain: {}'.format(i + 1, prot_df.loc[i, 'prot'], prot_df.loc[i, 'chain']))
    threedi_dict = pp.fpstr.threedi(
        prot_name=prot_df.loc[i, 'prot'],
        prot_chain=prot_df.loc[i, 'chain'],
        pdb_fp=to('data/pdb/pdbtm/'),
        mode='chain',
    )
    print(threedi_dict)
```

:material-note-multiple-outline: Output
``` text
No.1: protein: 1aij chain: L
{'1aij': {'L': {'state': masked_array(data=[--, 2, 2, 12, 12, 5, 12, 17, 9, 9, 2, 0, 2, 9, 0, 12, 12, 2, 5, 1, 12, 13, 15, ..., 12, 14, --],
mask=[ True, False, False, False, False, False, False, False,
False, False, False, False, ..., False, False, False, True],
fill_value=2,
dtype=uint8), 'encoded_sequence': 'DDDPPGPVLLDADLAPPDGCPQSQADVPHRAGPLNVQLVVLVVVLVVLLVVLCVVVPHDDQQPRKQAAAPCVQFPHADDSNRHNSVVSSVVSVLRNQLSVLVSLSSVCSSVVHDSLPSVLSVLVSVLLCLQQPLLSVLLRHRRSHFMDGDVRSVVSVVVLCPQQPNVCLQVLLVLLVVLVVVLVVLVVLVVVLQVCQQPPPPPDDRHDQVVSQVVCCVVPNGDCGPVRSVVSNSVSNSSSSVSSSVSSSCDNNVDRHRSVVVVCVVCCPPPNVPDDDDPRD'}}}
No.2: protein: 1aig chain: L
{'1aig': {'L': {'state': masked_array(data=[--, 2, 2, 12, ..., 2, 2, 11, --],
mask=[ True, False, False, ..., False, False, False,
True],
fill_value=2,
dtype=uint8), 'encoded_sequence': 'DDDPPCVVLLDADQADPDHCPQQDDDVPHGQGPLNVQLVVLVVVLVVLLVVQCVVVPHDDQQPSKQAAADPVCWADADDPNNHRSSVSSVVSVLSNQLSVLVVVSSVCRSVVHDSPVSVLSVLVSVLQCLQQPVLSVLLGHNRSHFMPGDVRSVVSVVVQCPLQPNVCLQVLLVVLVVLVVVLVVLVVVVVVLQCCQQVPPPPDDRHDQVVSQVVCCVVPVGDCGPVSSVVSNSCSSSSSSVSNSVSVSCEPRVDRHRSVVVCCCVQCVPPNNPDDDDDND'}}}
No.3: protein: 1xqf chain: A
{'1xqf': {'A': {'state': masked_array(data=[--, 0, 0, ..., 14, 7,
5, 2, --],
mask=[ True, False, False, ..., False,
False,  True],
fill_value=2,
dtype=uint8), 'encoded_sequence': 'DAADVQLQVLLVVLLVLLLLLQVPQLLLQLLQQDDVVQNVQLVVLLVVVLVVLLVVLQQAQVQFQAACDAQFTHDRPQGNHPPQDQRDDDPRHGVVSVSSNVSSLLSVLLSLQSSVCSVWWDSVLSVLLSVLCSVQALRRLSCNCPRPHVLVVLQAFFQQCLLNRQQLSLLLNVVVVVCHNHLVSSVVSLVSNLRSLLSRQLCSVSGPDVLSVQLNLLLVQLLVLQLVLQQVVCCVPVVGGDSVSSSQSSLQSSSLCRRARSWAGSVLSNVSSNVSSNQLNVQLVVVSSVSSSSNRSSLSSLQSNLPGSPVVNGTVHGPPPDDSVSNNVSSVVSSVVSNVSSNVSSVVSSVVSCVPPRIGDD'}}}
```