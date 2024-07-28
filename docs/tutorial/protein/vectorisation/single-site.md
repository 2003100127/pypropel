1. install TMKit tool,

Read a protein sequence.

:material-language-python: Python
``` py linenums="1"
import tmkit as tmk

sequence = tmk.seq.read_from_fasta(
    fasta_fpn='./data/fasta/1xqfA.fasta'
)
print(sequence)
```

:material-note-multiple-outline: Output
``` shell
ALLSFERKYRVPGGTLVGGNLFDFWVGPFYVGFFGVATFFFAALGIILIAWSAVLQGTWNPQLISVYPPALEYGLGGAPLAKGGLWQIITICATGAFVSWALREVEICRKLGIGYHIPFAFAFAILAYLTLVLFRPVMMGAWGYAFPYGIWTHLDWVSNTGYTYGNFHYNPAHMIAISFFFTNALALALHGALVLSAANPEKGKEMRTPDHEDTFFRDLVGYSIGTLGIHRLGLLLSLSAVFFSALCMIITGTIWFDQWVDWWQWWVKLPWWANIPGGING
```

:material-language-python: Python
``` py linenums="1"
pos_list = tmk.seq.pos_list_single(
    len_seq=len(sequence),
    seq_sep_superior=None,
    seq_sep_inferior=0,
)
print(pos_list)
```

:material-note-multiple-outline: Output
``` shell
[
    [1],
    [2],
    [3],
    ...,
    [280],
    [281],
]
```


:material-language-python: Python
``` py linenums="1"
positions = tmk.seq.pos_single(sequence=sequence, pos_list=pos_list)
print(positions)
```

:material-note-multiple-outline: Output
``` shell
[
    [1, 'A', 1, 0],
    [2, 'L', 2, 0],
    [3, 'L', 3, 0],
    ...,
    [280, 'N', 280, 0],
    [281, 'G', 281, 0],
]
```


:material-language-python: Python
``` py linenums="1"
win_aa_ids = tmk.seq.win_id_single(
    sequence=sequence,
    position=positions,
    window_size=1,
)
print(win_aa_ids)
```

:material-note-multiple-outline: Output
``` shell
[
    [None, 1, 2],
    [1, 2, 3],
    [2, 3, 4],
    ..., 
    [279, 280, 281],
    [280, 281, None],
]
```

:material-language-python: Python
``` py linenums="1"
win_aas = tmk.seq.win_name_single(
    sequence=sequence,
    position=positions,
    window_size=1,
    mids=win_aa_ids,
)
print(win_aas)
```

:material-note-multiple-outline: Output
``` shell
[
    [None, 'A', 'L'],
    ['A', 'L', 'L'],
    ['L', 'L', 'S'],
    ...,
    ['G', 'I', 'N'],
    ['I', 'N', 'G'],
    ['N', 'G', None],
]
```

:material-language-python: Python
``` py linenums="1"
features = [[] for i in range(len(sequence))]
print(features)
print(len(features))
```

:material-note-multiple-outline: Output
``` shell
[[], [], [], ..., [], [], []]
281
```
