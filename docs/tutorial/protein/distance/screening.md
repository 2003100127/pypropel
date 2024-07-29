## 1. Label residues

The following code allows users to label residues as interation sites or non-interaction sites. The labels will be rendered as in column `is_contact`.

`dist_fp` will be the directory where the `*.dist` file is stored.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

pp.dist.labelling(
    dist_fp=to('data/pdb/complex/pdbtm/'),
    prot_name='1aij',
    file_chain='L',
    cutoff=6,
)
```

:material-note-multiple-outline: Output
``` shell
28/07/2024 21:33:30 logger: ================>Labeling data...
28/07/2024 21:33:30 logger: ================>Time to label distances 1aij L: 0.004000186920166016s.
     fasta_id aa  pdb_id    dist_1     dist_2  is_contact
0           1  A       1  3.359507   3.820152           1
1           2  L       2  4.832898   2.840869           1
2           3  L       3  3.621853   3.450767           1
3           4  S       4  6.046180   2.841488           1
4           5  F       5  3.438960   3.688982           1
..        ... ..     ...       ...        ...         ...
276       277  G     277  2.910155  43.177975           1
277       278  G     278  3.036153  40.448544           1
278       279  I     279  2.979592  36.381110           1
279       280  N     280  2.980944  38.922234           1
280       281  G     281  3.484761  42.785560           1

[281 rows x 6 columns]
```

## 2. Find interaction partners

The following code allows you to find the partner chains in interaction with the target protein chain within one protein complex.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

pp.dist.interation_partners(
    dist_fp=to('data/pdb/complex/pdbtm/'),
    prot_name='1aij',
    file_chain='L',
    cutoff=6,
    pdb_fp=to('data/pdb/complex/pdbtm/'),
)
```

It will output the following list containing the partner chains interacting with chain `L` of protein `1aij`.

:material-note-multiple-outline: Output
``` shell
['M', 'H']
```

