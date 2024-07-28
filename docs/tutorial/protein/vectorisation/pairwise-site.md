1. install TMKit tool,

Read a protein sequence.

:material-language-python: Python
``` py linenums="1"
import tmkit as tmk

sequence = sfasta().get(
    fasta_fpn=to("data/fasta/1aigL.fasta")
)
print(sequence)
```

:material-note-multiple-outline: Output
``` shell
ALLSFERKYRVPGGTLVGGNLFDFWVGPFYVGFFGVATFFFAALGIILIAWSAVLQGTWNPQLISVYPPALEYGLGGAPLAKGGLWQIITICATGAFVSWALREVEICRKLGIGYHIPFAFAFAILAYLTLVLFRPVMMGAWGYAFPYGIWTHLDWVSNTGYTYGNFHYNPAHMIAISFFFTNALALALHGALVLSAANPEKGKEMRTPDHEDTFFRDLVGYSIGTLGIHRLGLLLSLSAVFFSALCMIITGTIWFDQWVDWWQWWVKLPWWANIPGGING
```

:material-language-python: Python
``` py linenums="1"
pos_list = tmk.seq.pos_list_pair(
    len_seq=len(sequence),
    seq_sep_superior=None,
    seq_sep_inferior=0,
)
print(pos_list)
```

:material-note-multiple-outline: Output
``` shell
[
    [1, 2],
    [1, 3],
    [1, 4],
    ...,
    [279, 281],
    [280, 281],
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
    [1, 'A', 1, 2, 'L', 2, 0],
    [1, 'A', 1, 3, 'L', 3, 0],
    [1, 'A', 1, 4, 'S', 4, 0],
    ...,
    [279, 'I', 279, 281, 'G', 281, 0],
    [280, 'N', 280, 281, 'G', 281, 0],
]
```


:material-language-python: Python
``` py linenums="1"
window_size = 1
win_aa_ids = tmk.seq.win_id_pair(
    sequence=sequence,
    position=positions,
    window_size=window_size,
)
print(win_aa_ids[:3])
```

:material-note-multiple-outline: Output
``` shell
[
    [[None, 1, 2], [1, 2, 3]],
    [[None, 1, 2], [2, 3, 4]],
    [[None, 1, 2], [3, 4, 5]],
    ...,
    [[278, 279, 280], [280, 281, None]],
    [[279, 280, 281], [280, 281, None]],
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
    [None, 'A', 'A', 'L', 'L', 'L'],
    [None, 'L', 'A', 'L', 'L', 'S'],
    [None, 'L', 'A', 'S', 'L', 'F'],
    ...,
    ['G', 'N', 'I', 'G', 'N', None],
    ['I', 'N', 'N', 'G', 'G', None],
]
```

:material-language-python: Python
``` py linenums="1"
features = [[] for i in range(len(win_aa_ids))]
print(features)
print(len(features))
```

:material-note-multiple-outline: Output
``` shell
[[], [], [], [], ..., [], [], []]
39340
```
