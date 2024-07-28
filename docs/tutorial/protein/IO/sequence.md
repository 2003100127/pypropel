## Read
PyPropel reads a protein sequence in FASTA format with the `pp.seq.read` module.

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

## Save
A protein sequence can be saved.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

seq = pp.seq.save(
    list_2d=[
        ['1aigL-new', seq],
    ],
    sv_fp=to('data/fasta/'),
)
```

:material-note-multiple-outline: Output
``` shell
No.1 saving 1aigL-new in FASTA format.
```

## Check empty
You can check if all sequences are empty.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

prot_df = pd.DataFrame({
    'prot': ['1aig', '1aij', '1xqf'],
    'chain': ['L', 'L', 'A'],
})

seq = pp.seq.is_empty(
    prot_df=prot_df,
    fasta_fp=to('data/fasta/'),
    sv_empty_fp=to('data/'),
)
```

:material-note-multiple-outline: Output
``` shell
27/07/2024 22:14:29 logger: ============>No0. protein 1aig chain L
27/07/2024 22:14:29 logger: ============>No1. protein 1aij chain L
27/07/2024 22:14:29 logger: ============>No2. protein 1xqf chain A
Finished
```

## Check match
You can check if all sequences in one format match those in another format. Below shows if sequences in FASTA format match those in XML format.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

prot_df = pd.DataFrame({
    'prot': ['1aig', '1aij', '1xqf'],
    'chain': ['L', 'L', 'A'],
})

seq = pp.seq.is_match(
    prot_df=prot_df,
    fasta_path=to('data/fasta/'),
    # pdb_path=to('data/pdb/pdbtm/'),
    xml_path=to('data/xml/'),
    kind='fasta<->xml',
    sv_mismatch_fp=to('data/'),
)
```

:material-note-multiple-outline: Output
``` shell
27/07/2024 22:17:32 logger: ============>No0. protein 1aig chain L
27/07/2024 22:17:32 logger: ============>They match.
27/07/2024 22:17:32 logger: ============>No1. protein 1aij chain L
27/07/2024 22:17:32 logger: ============>They match.
27/07/2024 22:17:32 logger: ============>No2. protein 1xqf chain A
27/07/2024 22:17:32 logger: ============>They do not match.
Finished
```