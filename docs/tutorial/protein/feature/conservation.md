## Entropy

In the context of multiple sequence alignments (MSA), entropy measures the variability or uncertainty at each position in the alignment. It quantifies how conserved or variable a particular position is across the aligned sequences, offering insight into which positions are likely to be functionally or structurally important.

We can first read a MSA in ALN format with `pp.msa.msaparser`. Then, taking this MSA as input to `pp.fpmsa.entropy`, we obtain a dictionary of entropy values.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

msa = pp.msa.msaparser(msa_fpn=to('data/msa/aln/1aijL.aln')).read()

entropy_dict = pp.fpmsa.entropy(msa=msa)

print(entropy_dict)
```

:material-note-multiple-outline: Output
``` shell
{1: 0.1960388659924533, 2: 0.22426701608683547, 3: 0.1850087818623904, 4: 0.20286277711377826, 5: 0.18773405896608028, 6: 0.18177904812214024, 7: 0.22124858188044072, 8: 0.1964272886606503, 9: 0.18956787701644068, 10: 0.19008797269348587, ..., 279: 0.07534724091084649, 280: 0.07102837084165835, 281: 0.06073066166727346}
```

If we want to consider the information from a gap symbol '-', we can achieve it with `pp.fpmsa.entropy_gap`.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

entropy_dict = pp.fpmsa.entropy_gap(msa=msa, gap_thres=100)

print(entropy_dict)
```

:material-note-multiple-outline: Output
``` shell
{1: 0.3256135083150946, 2: 0.3724994507477074, 3: 0.30729293513478106, ..., 279: 0.12514905822699368, 280: 0.11797557031649158, 281: 0.10087144560680732}
```

## Custom-conservation

We borrow a technical idea from [^1] to build a custom-conservation score based on the above calculated entropy values.

[^1]: Zeng B, Hönigschmid P, Frishman D. Residue co-evolution helps predict interaction sites in α-helical membrane proteins. J Struct Biol. 2019 May 1;206(2):156-169. doi: 10.1016/j.jsb.2019.02.009. Epub 2019 Mar 2. PMID: 30836197. https://doi.org/10.1016/j.jsb.2019.02.009

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

cc = pp.fpmsa.conservation_custom(ent_dict=entropy_dict)

print(cc)
```

:material-note-multiple-outline: Output
``` shell
{1: 0.891307540667243, 2: 0.8756566285859075, 3: 0.8974230982362708, ..., 279: 0.9582242180545316, 280: 0.9606187871466461, 281: 0.9663282842404544}
```

## Mutual information

Mutual information (MI) is an information theory metric that quantifies how much information one random variable reveals about another. In the context of multiple sequence alignments (MSA), MI measures the statistical dependency between sequence positions, helping to identify co-evolving residues. 

To calculate the MI between columns 1 and 2, we can do

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

msa = pp.msa.msaparser(msa_fpn=to('data/msa/aln/1aijL.aln')).read()

mi = pp.fpmsa.mutual_information(msa=msa, i=1, j=2)

print(mi)
```

:material-note-multiple-outline: Output
``` shell
0.20924659742293167
```


## Jensen–Shannon Divergence

Jensen–Shannon Divergence (JSD) [^3] is a useful metric in structual bioinformatics for evaluating the conservation of sequence positions in multiple sequence alignments (MSAs). By analyzing the distribution of amino acids or nucleotides at each position, JSD helps to identify conserved and variable regions, indicating areas of functional or structural significance.

[^3]: John A. Capra, Mona Singh, Predicting functionally important residues from sequence conservation, Bioinformatics, Volume 23, Issue 15, August 2007, Pages 1875–1882, https://doi.org/10.1093/bioinformatics/btm270

A JSD file looks like
``` shell
# ./CLEC2B_LOC113845378.clustal -- js_divergence - window_size: 3 - window lambda: 0.50 - background: blosum62 - seq. weighting: True - gap penalty: 1 - normalized: False
        # align_column_number	score	column
        0	-1000.000000	M-M-M-T----TET--TTTMM-M-----------------------------------------------------P-----------------------------------------------
        1	-1000.000000	E-S-D-Q----QTE--QQQEE-E-----------------------------------------------------V----------------------G------------------------
        2	-1000.000000	K-S-S-N----DGY--DNNPP-P-----------------------------------------------------P----------------------P------------------------
        3	-1000.000000
        ...
        205	-1000.000000	PN---S---P--VPF------R---------------L---LL---LLL------L----DD---V----A--E----------------------------------EME---L-EE-I----
        206	-1000.000000	EP---P---R---VP----------------------S---SS---SSS------S----HH---T----E--K----------------------------------SCS---C-SS-M----
```

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

jsd = pp.fpmsa.jsd(
    fpn=to('data/conservation/jsd/SR24_CtoU/CLEC2B_LOC113845378.jsd'),
    mode='standalone',
)
print(jsd)
```

The `standalone` mode means the output returned from running its standalone package.

:material-note-multiple-outline: Output
``` shell
     alignment_col  score                                                seq
0                0    0.0  M-M-M-T----TET--TTTMM-M-----------------------...
1                1    0.0  E-S-D-Q----QTE--QQQEE-E-----------------------...
2                2    0.0  K-S-S-N----DGY--DNNPP-P-----------------------...
3                3    0.0  E-E-E-E----EEN--EEEAA-A-----------------------...
4                4    0.0  V-N-N-G----EKN--EGG---------------------------...
..             ...    ...                                                ...
202            202    0.0  KL---T---G--TTT------E-P---------D--DEGDDEES--...
203            203    0.0  II---P---S--PPP------Y-L---------S--SQ-QSRRK--...
204            204    0.0  SV---F---V--PVF------M-L-------------F-S-LL---...
205            205    0.0  PN---S---P--VPF------R---------------L---LL---...
206            206    0.0  EP---P---R---VP----------------------S---SS---...

[207 rows x 3 columns]
```


## ConSurf

ConSurf [^2] is designed to identify and visualize conserved regions within protein or nucleic acid sequences. By mapping conservation levels onto the three-dimensional structure of a protein, ConSurf aids in pinpointing functionally and structurally significant regions.

[^2]: Barak Yariv, Elon Yariv, Amit Kessel, Gal Masrati, Adi Ben Chorin, Eric Martz, Itay Mayrose, Tal Pupko, and Nir Ben-Tal
Using evolutionary data to make sense of macromolecules with a 'face-lifted' ConSurf. Protein Science 2023; DOI: 10.1002/pro.4582; PMID: 36718848. https://doi.org/10.1002/pro.4582

We tease out the information from it by

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

consurf = pp.fpmsa.consurf(
    fpn=to('data/conservation/consurf/E.consurf'),
    mode='v1'
)
print(consurf)
```

The `v1` mode means the output adopted by ConSurf before 2024.

:material-note-multiple-outline: Output
``` shell
    position amino acid  score color exposed/buried structral/functional
0          1          M -1.409     9              e                    f
1          2          Y  0.136    5*              e                     
2          3          S -0.426     6              e                     
3          4          F -0.115     5              b                     
4          5          V -1.201     8              b                     
..       ...        ...    ...   ...            ...                  ...
70        71          P  0.637    3*              e                     
71        72          D  0.109    5*              e                     
72        73          L -0.200    6*              b                     
73        74          L -0.120    5*              e                     
74        75          V -0.991     8              e                    f

[75 rows x 6 columns]
```

The original output of ConSurf of [E protein](https://www.ncbi.nlm.nih.gov/gene/43740570) is
``` shell
    0  1       2  3   4              5  6   ...   8  9  10 11 12      13         14
0    1  M  -1.409      9  -1.875,-1.232     ...  9,8        e  f  78/150        M,N
1    2  Y   0.136     5*  -0.618, 0.577     ...  7,3        e     78/150    C,F,Y,L
2    3  S  -0.426      6  -1.003,-0.074     ...  8,5        e     78/150  E,P,S,Y,D
3    4  F  -0.115      5  -0.776, 0.279     ...  7,4        b     78/150      I,F,L
4    5  V  -1.201      8  -1.582,-0.928     ...  9,7        b     78/150      V,F,Q
..  .. ..     ... ..  ..            ... ..  ...  ... .. .. .. ..     ...        ...
70  71  P   0.637     3*  -0.365, 1.225     ...  6,2        e     34/150      S,P,L
71  72  D   0.109     5*  -0.776, 0.755     ...  7,3        e     34/150        D,E
72  73  L  -0.200     6*  -1.003, 0.419     ...  8,4        b     34/150        L,F
73  74  L  -0.120     5*  -1.003, 0.419     ...  8,4        e     34/150        I,L
74  75  V  -0.991      8  -1.582,-0.618     ...  9,7        e  f  33/150          V

[75 rows x 15 columns]
```

The columns from the left to the right are defined as:

!!! info "Column definition" 

    ``` shell
    - POS: The position of the AA in the SEQRES derived sequence.
    - SEQ: The SEQRES derived sequence in one letter code.
    - SCORE: The normalized conservation scores.
    - COLOR: The color scale representing the conservation scores (9 - conserved, 1 - variable).
    - CONFIDENCE INTERVAL: When using the bayesian method for calculating rates, a confidence interval is assigned to each of the inferred evolutionary conservation scores.
    - CONFIDENCE INTERVAL COLORS: When using the bayesian method for calculating rates. The color scale representing the lower and upper bounds of the confidence interval.
    - B/E: Burried (b) or Exposed (e) residue.
    - FUNCTION: functional (f) or structural (s) residue (f - highly conserved and exposed, s - highly conserved and burried).
    - MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position.
    - RESIDUE VARIETY: The residues variety at each position of the multiple sequence alignment.
    
    ```