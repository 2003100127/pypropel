## Absolute position

We can use PyPropel to calculate the positional features of a residue in a protein.

First, PyPropel needs to read a protein sequence like below.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

seq = pp.seq.read(
    fasta_fpn=to('data/fasta/1aigL.fasta'),
)
print(seq)
```

:material-note-multiple-outline: Output
``` shell
ALLSFERKYRVPGGTLVGGNLFDFWVGPFYVGFFGVATFFFAALGIILIAWSAVLQGTWNPQLISVYPPALEYGLGGAPLAKGGLWQIITICATGAFVSWALREVEICRKLGIGYHIPFAFAFAILAYLTLVLFRPVMMGAWGYAFPYGIWTHLDWVSNTGYTYGNFHYNPAHMIAISFFFTNALALALHGALVLSAANPEKGKEMRTPDHEDTFFRDLVGYSIGTLGIHRLGLLLSLSAVFFSALCMIITGTIWFDQWVDWWQWWVKLPWWANIPGGING
```

Then, we use `pp.fpsite.pos_abs_val` to calculate a residue positioned at 1 in the protein.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

p = pp.fpsite.pos_abs_val(
    pos=1,
    seq=sequence,
)
print(p)
```

:material-note-multiple-outline: Output
``` shell
0.0035587189
```

## Relative position

We implement a function to calculate the positional feature of a residue relative to the transmembrane segment, proposed by Zeng et al [^1], which is to calculate the position of a residue relative to the boundary of a transmembrane segment (interval).

[^1]: Zeng B, Hönigschmid P, Frishman D. Residue co-evolution helps predict interaction sites in α-helical membrane proteins. J Struct Biol. 2019 May 1;206(2):156-169. doi: 10.1016/j.jsb.2019.02.009. Epub 2019 Mar 2. PMID: 30836197. https://doi.org/10.1016/j.jsb.2019.02.009

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

pos = pp.fpsite.pos_rel_val(
    pos=5,
    interval=[4, 10],
)
print(pos)
```

:material-note-multiple-outline: Output
``` shell
0.16666666666666666
```

## One-hot position

In addition, positions of residues can be one-hot encoded to describe their rough positional occurrence in the protein. There are two residue contact prediction tools that use this idea.

#### MetaPSICOV

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

pos = pp.fpsite.metapsicov()
print(pos)
```

:material-note-multiple-outline: Output
``` shell
{
    'v<5': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
    'v=5': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
    'v=6': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    'v=7': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    'v=8': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    'v=9': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'v=10': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'v=11': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'v=12': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'v=13': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    '14<=v<18': [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    '18<=v<23': [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    '23<=v<28': [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    '28<=v<38': [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    '38<=v<48': [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'v>=48': [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'none': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
}
```


#### DeepConPred

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

pos = pp.fpsite.deepconpred()
print(pos)
```

:material-note-multiple-outline: Output
``` shell
{
    '24<=v<=28': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
 
    '29<=v<=33': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0],
    '34<=v<=38': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0],
    '39<=v<=43': [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    '44<=v<=48': [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    '49<=v<=58': [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    '59<=v<=68': [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    '69<=v<=78': [0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    '79<=v<=88': [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'v>=89': [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    'none': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
}
```