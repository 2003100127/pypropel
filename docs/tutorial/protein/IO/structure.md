## 1. Read
PyPropel pre-processes a protein structure in PDB format with a variety of functions to generate a new structure tailored for the analysis need.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

seq = pp.str.read(
    pdb_path=to('data/pdb/pdbtm/'),
    pdb_name='1aij',
    file_chain='L',
    seq_chain='L',
)
print(seq)
```

:material-note-multiple-outline: Output
``` shell
ALLSFERKYRVPGGTLVGGNLFDFWVGPFYVGFFGVATFFFAALGIILIAWSAVLQGTWNPQLISVYPPALEYGLGGAPLAKGGLWQIITICATGAFVSWALREVEICRKLGIGYHIPFAFAFAILAYLTLVLFRPVMMGAWGYAFPYGIWTHLDWVSNTGYTYGNFHYNPAHMIAISFFFTNALALALHGALVLSAANPEKGKEMRTPDHEDTFFRDLVGYSIGTLGIHRLGLLLSLSAVFFSALCMIITGTIWFDQWVDWWQWWVKLPWWANIPGGING
```


## 2. Splitting a complex into chains 
If there is a list of protein chains, the following code allows you to extract single chains from a protein complex.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

prot_df = pd.DataFrame({
    'prot': ['1aig', '1aij', '1xqf'],
    'chain': ['L', 'L', 'A'],
})
    
pp.str.split_cplx_to_sgl(
    prot_df=prot_df,
    pdb_path=to('data/pdb/complex/pdbtm/'),
    sv_fp=to('data/'),
)
```

:material-note-multiple-outline: Output
``` shell
28/07/2024 17:28:16 logger: ============>No0. protein 1aig chain L
D:\Programming\anaconda3\envs\prot\Lib\site-packages\Bio\PDB\Atom.py:232: PDBConstructionWarning: Could not assign element 'M' for Atom (name=MG) with given element 'M'
  warnings.warn(msg, PDBConstructionWarning)
28/07/2024 17:28:16 logger: ================>success in building 1aigL model.
28/07/2024 17:28:16 logger: ============>No1. protein 1aij chain L
D:\Programming\anaconda3\envs\prot\Lib\site-packages\Bio\PDB\Atom.py:232: PDBConstructionWarning: Could not assign element 'M' for Atom (name=MG) with given element 'M'
  warnings.warn(msg, PDBConstructionWarning)
28/07/2024 17:28:16 logger: ================>success in building 1aijL model.
28/07/2024 17:28:16 logger: ============>No2. protein 1xqf chain A
28/07/2024 17:28:16 logger: ================>success in building 1xqfA model.
Finished
```


## 3. Delete END from PDB
In reality, we might extract the 3D coordinates of a single protein chain from a PDB structure. Due to unproper formatting, the extracted structure may contain something irrelevant, posing a challenge for downstream analysis. For example, to get the correct information about the relative
solvent accessibility (RSA), DSSP needs a PDB file without `END` in the end of the PDB file. We can use the following code to remove the mark in bulk.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

prot_df = pd.DataFrame({
    'prot': ['1aig', '1aij', '1xqf'],
    'chain': ['L', 'L', 'A'],
})
    
pp.str.del_end(
    prot_df,
    sv_fp=to('data/'),
    pdb_path=to('data/pdb/pdbtm/'),
)
```


Before deletion, the PDB file has a `END` in the end.
``` shell
ATOM      1  N   ALA L   1     -27.710  -2.809 -20.227  1.00 61.86           N  
ATOM      2  CA  ALA L   1     -26.559  -3.484 -19.659  1.00 62.35           C  
ATOM      3  C   ALA L   1     -25.256  -3.284 -20.413  1.00 62.67           C  
......
ATOM   2232  OXT GLY L 281      28.878   2.108  26.859  1.00 74.20           O  
TER    2233      GLY L 281                                                       
END
```

:material-note-multiple-outline: Output
```shell
28/07/2024 17:03:30 logger: ============>No0. protein 1aig chain L
28/07/2024 17:03:30 logger: ===============>Successfully reformatted
28/07/2024 17:03:30 logger: ============>No1. protein 1aij chain L
28/07/2024 17:03:30 logger: ===============>Successfully reformatted
28/07/2024 17:03:30 logger: ============>No2. protein 1xqf chain A
28/07/2024 17:03:30 logger: ===============>Successfully reformatted
Finished
```

After `pp.str.del_end` is used, the PDB file looks like.
``` shell
ATOM      1  N   ALA L   1     -27.710  -2.809 -20.227  1.00 61.86           N  
ATOM      2  CA  ALA L   1     -26.559  -3.484 -19.659  1.00 62.35           C  
ATOM      3  C   ALA L   1     -25.256  -3.284 -20.413  1.00 62.67           C  
......
ATOM   2232  OXT GLY L 281      28.878   2.108  26.859  1.00 74.20           O  
TER    2233      GLY L 281                                                 
```


## 4. Remove HETATM

In molecular structures of PDB (Protein Data Bank) files, HETATM (Hetero Atom) refers to atoms that are not part of standard amino acid or nucleic acid residues. These heteroatoms typically include:

1. Ligands: Small molecules or ions that bind to proteins or nucleic acids.
2. Cofactors: Non-protein chemical compounds or metallic ions required for protein activity.
3. Water Molecules: Included in crystal structures.

To remove this part, users can do

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

prot_df = pd.DataFrame({
    'prot': ['1aig', '1aij', '1xqf'],
    'chain': ['L', 'L', 'A'],
})
    
pp.str.remove_hetatm(
    prot_df=prot_df,
    pdb_path=to('data/pdb/complex/pdbtm/'),
    sv_fp=to('data/pdb/'),
)
```

:material-note-multiple-outline: Output
```shell
28/07/2024 17:44:07 logger: ============>No.1 protein 1aig chain L
28/07/2024 17:44:07 logger: ============>No.2 protein 1aij chain L
28/07/2024 17:44:07 logger: ============>No.3 protein 1xqf chain A
```