## Relative solvent accessibility

Secondary structure prediction methods are used to predict the local structures that form within a protein, specifically alpha helices, beta strands, and loops or turns. These predictions are essential for understanding protein function and guiding experimental structure determination.

Secondary Structures:

1. Alpha Helices: Coiled structures stabilized by hydrogen bonds between the backbone atoms.
2. Beta Strands: Extended strands that can form hydrogen bonds with neighboring strands to create a sheet-like arrangement.
3. Loops/Turns: Irregular, flexible regions that connect helices and strands.

PyPropel can handle the output of 3 secondary structure prediction methods, `PSIPRED`[^1], `SSpro`[^2], and `Spider3`[^3].

[^1]: Jones DT. Protein secondary structure prediction based on position-specific scoring matrices. J Mol Biol. 1999 Sep 17;292(2):195-202. doi: 10.1006/jmbi.1999.3091. PMID: 10493868.

[^2]: Urban G, Magnan CN, Baldi P. SSpro/ACCpro 6: almost perfect prediction of protein secondary structure and relative solvent accessibility using profiles, deep learning and structural similarity. Bioinformatics. 2022 Mar 28;38(7):2064-2065. doi: 10.1093/bioinformatics/btac019. PMID: 35108364.

[^3]: Jack Hanson, Kuldip Paliwal, Thomas Litfin, Yuedong Yang, Yaoqi Zhou, Improving prediction of protein secondary structure, backbone angles, solvent accessibility and contact numbers by using predicted contact maps and an ensemble of recurrent and residual convolutional neural networks, Bioinformatics, Volume 35, Issue 14, July 2019, Pages 2403–2410, https://doi.org/10.1093/bioinformatics/bty1006

#### PSIPRED

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

df = pp.fpsite.ss_psipred(
    psipred_path=to('data/ss/psipred/'),
    prot_name='1aig',
    file_chain='L',
    kind='ss', # horiz, ss, ss2
)
print(df)
```

:material-note-multiple-outline: Output
``` shell
       0  1  2      3      4      5
0      1  A  C  0.992  0.003  0.004
1      2  L  C  0.697  0.082  0.148
2      3  L  C  0.684  0.150  0.112
3      4  S  C  0.673  0.121  0.155
4      5  F  C  0.363  0.245  0.198
..   ... .. ..    ...    ...    ...
276  277  G  C  0.806  0.137  0.082
277  278  G  C  0.765  0.123  0.104
278  279  I  C  0.733  0.113  0.133
279  280  N  C  0.825  0.088  0.068
280  281  G  C  0.997  0.003  0.002

[281 rows x 6 columns]
```

#### SSpro

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

df = pp.fpsite.ss_sspro(
    sspro_path=to('data/ss/sspro/'),
    prot_name='1aig',
    file_chain='L'
)
print(df)
```

:material-note-multiple-outline: Output
``` shell
     0
0    C
1    E
2    C
3    C
4    C
..  ..
276  C
277  C
278  C
279  C
280  C
```

#### SSpro8

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
0    C
1    B
2    C
3    T
4    T
..  ..
276  S
277  S
278  S
279  C
280  C

[281 rows x 1 columns]
```


#### SSpro8

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
0     C
1     C
2     C
3     C
4     C
     ..
70    C
71    C
72    C
73    C
74    C
Name: SS, Length: 75, dtype: object
```

The original output of Spider3 contains multiple columns to describe protein properties. `P(C)`, `P(H)`, and `P(E)` are three columns are predicted probabilities of three types of secondary structures `C`, `H`, and `E`.

``` text
     # SEQ SS    ASA    Phi  ...  HSE_beta_down     CN   P(C)   P(H)   P(E)
0    1   M  C  1.578 -0.945  ...          0.046  0.085  0.980  0.019  0.001
1    2   Y  C  1.640 -0.993  ...          0.038  0.098  0.931  0.060  0.009
2    3   S  C  0.873 -0.817  ...          0.102  0.144  0.868  0.127  0.005
3    4   F  C  1.542 -0.857  ...          0.068  0.120  0.744  0.253  0.004
4    5   V  C  1.032 -0.892  ...          0.070  0.124  0.620  0.377  0.003
..  ..  .. ..    ...    ...  ...            ...    ...    ...    ...    ...
70  71   P  C  0.940 -0.679  ...          0.099  0.159  0.843  0.148  0.008
71  72   D  C  1.142 -0.741  ...          0.113  0.157  0.726  0.243  0.031
72  73   L  C  1.205 -0.820  ...          0.060  0.126  0.708  0.245  0.047
73  74   L  C  1.180 -0.859  ...          0.075  0.139  0.824  0.109  0.068
74  75   V  C  1.153 -0.910  ...          0.075  0.150  0.983  0.006  0.011

[75 rows x 16 columns]
```