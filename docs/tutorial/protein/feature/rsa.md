## Relative solvent accessibility

Relative solvent accessibility (RSA) is a metric used in structural biology to quantify the extent to which an amino acid residue in a protein is exposed to solvent, typically water. RSA offers insights into the environmental context of the residue and its potential functional roles within the protein structure.

PyPropel can access the RSA of 2 methods, `solvpred`[^1] and `accpro`[^2].

[^1]: Jones DT. Protein secondary structure prediction based on position-specific scoring matrices. J Mol Biol. 1999 Sep 17;292(2):195-202. doi: 10.1006/jmbi.1999.3091. PMID: 10493868.

[^2]: Urban G, Magnan CN, Baldi P. SSpro/ACCpro 6: almost perfect prediction of protein secondary structure and relative solvent accessibility using profiles, deep learning and structural similarity. Bioinformatics. 2022 Mar 28;38(7):2064-2065. doi: 10.1093/bioinformatics/btac019. PMID: 35108364.

#### solvpred

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

df = pp.fpsite.rsa_solvpred(
    solvpred_fp=to('data/accessibility/solvpred/'),
    prot_name='1aig',
    file_chain='L',
)
print(df)
```

:material-note-multiple-outline: Output
``` shell
       0  1      2
0      1  A  0.855
1      2  L  0.458
2      3  L  0.308
3      4  S  0.420
4      5  F  0.219
..   ... ..    ...
276  277  G  0.397
277  278  G  0.336
278  279  I  0.215
279  280  N  0.548
280  281  G  0.912

[281 rows x 3 columns]
```

#### accpro

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

df = pp.fpsite.rsa_accpro(
    accpro_fp=to('data/accessibility/accpro/'),
    prot_name='1aig',
    file_chain='L',
)
print(df)
```

:material-note-multiple-outline: Output
``` shell
     0
0    e
1    -
2    -
3    -
4    -
..  ..
276  e
277  e
278  e
279  e
280  e

[281 rows x 1 columns]
```

#### accpro20

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

df = pp.fpsite.rsa_accpro20(
    accpro20_fp=to('data/accessibility/accpro20/'),
    prot_name='1aig',
    file_chain='L',
)
print(df)
```


:material-note-multiple-outline: Output
``` shell
        0
0    0.30
1    0.00
2    0.00
3    0.05
4    0.00
..    ...
276  0.35
277  0.30
278  0.45
279  0.35
280  0.90

[281 rows x 1 columns]
```