
We will need to install [TMKit](https://github.com/2003100127/tmkit)[^1] to read a protein sequence and create single-site positions, placed with windows. This will initiate the vector of features, and be prepared for being fed by site-wise features.

[^1]: Jianfeng Sun, Arulsamy Kulandaisamy, Jinlong Ru, M Michael Gromiha, Adam P Cribbs, TMKit: a Python interface for computational analysis of transmembrane proteins, Briefings in Bioinformatics, Volume 24, Issue 5, September 2023, bbad288, https://doi.org/10.1093/bib/bbad288

After installation, we first read `1aigL.fasta`.

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

Generation of all posible residue pairs

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

Adding amino acid types and IDs to the positions of all posible residue pairs.

:material-language-python: Python
``` py linenums="1"
positions = tmk.seq.pos_pair(sequence=sequence, pos_list=pos_list)
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

Applying a sliding window to each residue pair.

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

Initiating feature vector.

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
