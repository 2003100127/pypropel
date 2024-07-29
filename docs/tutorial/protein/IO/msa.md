## Read
Users can operate multiple sequence alignments (MSAs) with the `pp.msa.read` module.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

msa = pp.msa.read(
    msa_fpn=to('data/msa/aln/1aijL.aln'),
)
print(msa[:3])

```

It returns a 1D list. Each element represents a homolog of the query protein sequence (the 1st protein). The first 3 sequences in the MSA of protein `1aijL` is given below.

:material-note-multiple-outline: Output
``` shell
['ALLSFERKYRVPGGTLVGGNLFDFWVGPFYVGFFGVATFFFAALGIILIAWSAVLQGTWNPQLISVYPPALEYGLGGAPLAKGGLWQIITICATGAFVSWALREVEICRKLGIGYHIPFAFAFAILAYLTLVLFRPVMMGAWGYAFPYGIWTHLDWVSNTGYTYGNFHYNPAHMIAISFFFTNALALALHGALVLSAANPEKGKEMRTPDHEDTFFRDLVGYSIGTLGIHRLGLLLSLSAVFFSALCMIITGTIWFDQWVDWWQWWVKLPWWANIPGGING', '-------------------------------------AVFCALMGTALIIWNTPLGPTWNLWQISVNPPDVKYGLGFAPLAEGGIWQWVTIFAIGAFCSWALREVEICRKLGIGYHVPFAFSFAIFAYVTLVVIRPVLMGSWSYGFPYGIFTHLDWVSNTGYQYGQFHWNPGHMIAITFFFTTCLALALHGGLVLSAINPDRGEPVKSPEHENTVFRDLVGYSIGTIGIHRVGLFLALSAVFWSAVCMLISGPVLGGSWPEWWEWWRRIPIWNP-------', '--------------------------------------IFFASLGICFIGYAASQGPTWDPFAISINPPDLKYGFAAAPLLEGGFWQAITVCALGAFFSWMLREVEISRKLGMGWHVPLAFCVPIFMFCVLQVFRPILMGGWGFAFPYGILSHLDWVNNFGFQYLNWHYNPGHMSSVSFLFCNAMALGLHGGLILSVTNPGDGDKVKTAEHENAYFRDVVGYSIGALAIHRLGLFLASNIFLTGAFGTIASGPFWTRGWPEWWGWWLDIPFWS--------',
...
]
```

## Frequency

We next calculate the frequencies of amino acids at each column.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

freq_col = pp.msa.freq_col_sgl(msa=msa)
print(freq_col)
```

:material-note-multiple-outline: Output
``` shell
{'A': array([4.80810828e-02, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, 3.09004388e-04, 0.00000000e+00,
       0.00000000e+00, 6.18008776e-05, ..., 1.85402633e-04, 6.79809653e-04,
       6.18008776e-05]), 'C': array([0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00,
       0.00000000e+00, 0.00000000e+00, ...
       }
```

The frequency of amino acids across the whole MSA.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

freq_whole = pp.msa.freq_whole_sgl(msa=msa)
print(freq_whole)
```

:material-note-multiple-outline: Output
``` shell
{'A': 0.05826085292688736, 'C': 0.004700385606685579, 'D': 0.008600878716107662, 'E': 0.015672130729309737, 'F': 0.05294113895278523, 'G': 0.06534156201388167, 'H': 0.019759785927038455, 'I': 0.04494023459261235, 'K': 0.0027040633087310126, 'L': 0.05968535215833517, 'M': 0.022540825417799226, 'N': 0.01639108826946766, 'P': 0.03155605592517563, 'Q': 0.011760641022454832, 'R': 0.017077715813173088, 'S': 0.03935616241622517, 'T': 0.02184298134471232, 'V': 0.040838283818220966, 'W': 0.01871466930702302, 'Y': 0.0222107075628659, '-': 0.42443809036607894}
```

Then, we can caluclate the frequency of amino acid pairs for any two columns in the MSA.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

freq_whole = pp.msa.freq_whole_sgl(msa=msa)
print(freq_whole)
```

:material-note-multiple-outline: Output
``` shell
{'AA': 0.0, 'AC': 0.0, 'AD': 0.0, ..., '-W': 0.0, '-Y': 0.0, '--': 0.9465422408998208}
```

## Count

The count of amino acids at each column.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

cnt_col = pp.msa.cnt_col_sgl(msa=msa)
print(cnt_col)
```

:material-note-multiple-outline: Output
``` shell
{'A': [778, 0, 0, 0, 0, 0, 5, ...], 'C': [0, 0, 0, 0, 0, 0, 0, ...], ...}
```

## Representation

#### Onehot

We use one-hot encodings to represent a MSA matrix.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

oh_rep = pp.msa.representation_onehot(
    msa=msa # 1aijL
)
print(oh_rep)
print(oh_rep.shape)
```

The length of the column is 281*21=5901, where 21 means 21 symbols including gap `-`.

:material-note-multiple-outline: Output
``` text
[[1 0 0 ... 0 0 0]
 [0 0 0 ... 0 0 1]
 [0 0 0 ... 0 0 1]
 ...
 [0 0 0 ... 0 0 1]
 [0 0 0 ... 0 0 1]
 [0 0 0 ... 0 0 1]]
(16181, 5901)
```

#### Frequency

We use frequencies of amino acids to represent a MSA matrix.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

oh_rep = pp.msa.representation_freq(
    msa=msa # 1aijL
)
print(oh_rep)
print(oh_rep.shape)
```

:material-note-multiple-outline: Output
``` text
[[0.0481 0.0256 0.0526 ... 0.0024 0.0067 0.0085]
 [0.     0.     0.     ... 0.     0.     0.    ]
 [0.     0.     0.     ... 0.     0.     0.    ]
 ...
 [0.     0.     0.     ... 0.     0.     0.    ]
 [0.     0.     0.     ... 0.     0.     0.    ]
 [0.     0.     0.     ... 0.     0.     0.    ]]
 
```

## Split

We can split a string of a sequence into capital letters of a list.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

l = pp.msa.split(
    msa=msa # 1aijL
)
print(l)
```

:material-note-multiple-outline: Output
``` text
[['A', 'L', 'L', 'S', 'F', 'E', 'R', 'K', 'Y', 'R', 'V', 'P', 'G', 'G', 'T', 'L', 'V', 'G', 'G', 'N', 'L', 'F', 'D', 'F', 'W', 'V', 'G', 'P', 'F', 'Y', 'V', 'G', 'F', 'F', 'G', 'V', 'A', 'T', 'F', 'F', 'F', 'A', 'A', 'L', 'G', 'I', 'I', 'L', 'I', 'A', 'W', 'S', 'A', 'V', 'L', 'Q', 'G', 'T', 'W', 'N', 'P', 'Q', 'L', 'I', 'S', 'V', 'Y', 'P', 'P', 'A', 'L', 'E', 'Y', 'G', 'L', 'G', 'G', 'A', 'P', 'L', 'A', 'K', 'G', 'G', 'L', 'W', 'Q', 'I', 'I', 'T', 'I', 'C', 'A', 'T', 'G', 'A', 'F', 'V', 'S', 'W', 'A', 'L', 'R', 'E', 'V', 'E', 'I', 'C', 'R', 'K', 'L', 'G', 'I', 'G', 'Y', 'H', 'I', 'P', 'F', 'A', 'F', 'A', 'F', 'A', 'I', 'L', 'A', 'Y', 'L', 'T', 'L', 'V', 'L', 'F', 'R', 'P', 'V', 'M', 'M', 'G', 'A', 'W', 'G', 'Y', 'A', 'F', 'P', 'Y', 'G', 'I', 'W', 'T', 'H', 'L', 'D', 'W', 'V', 'S', 'N', 'T', 'G', 'Y', 'T', 'Y', 'G', 'N', 'F', 'H', 'Y', 'N', 'P', 'A', 'H', 'M', 'I', 'A', 'I', 'S', 'F', 'F', 'F', 'T', 'N', 'A', 'L', 'A', 'L', 'A', 'L', 'H', 'G', 'A', 'L', 'V', 'L', 'S', 'A', 'A', 'N', 'P', 'E', 'K', 'G', 'K', 'E', 'M', 'R', 'T', 'P', 'D', 'H', 'E', 'D', 'T', 'F', 'F', 'R', 'D', 'L', 'V', 'G', 'Y', 'S', 'I', 'G', 'T', 'L', 'G', 'I', 'H', 'R', 'L', 'G', 'L', 'L', 'L', 'S', 'L', 'S', 'A', 'V', 'F', 'F', 'S', 'A', 'L', 'C', 'M', 'I', 'I', 'T', 'G', 'T', 'I', 'W', 'F', 'D', 'Q', 'W', 'V', 'D', 'W', 'W', 'Q', 'W', 'W', 'V', 'K', 'L', 'P', 'W', 'W', 'A', 'N', 'I', 'P', 'G', 'G', 'I', 'N', 'G'], ...]
```