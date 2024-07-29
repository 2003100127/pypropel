## 1. Distance calculation

We calculate the distance between residues from every two chains within a protein complex. The `pp.dist.all_vs_all` makes a all-against-all comparison between residues.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

df = pp.dist.all_vs_all(
    pdb_fp=to('data/pdb/complex/pdbtm/'),
    pdb_name='1aij',
)
print(df)
```

:material-note-multiple-outline: Output
``` shell
D:\Programming\anaconda3\envs\prot\Lib\site-packages\Bio\PDB\Atom.py:232: PDBConstructionWarning: Could not assign element 'M' for Atom (name=MG) with given element 'M'
  warnings.warn(msg, PDBConstructionWarning)
28/07/2024 18:41:22 logger: ============>chain 1 L and chain2 M
D:\Programming\anaconda3\envs\prot\Lib\site-packages\Bio\PDB\Polypeptide.py:144: BiopythonDeprecationWarning: 'three_to_one' will be deprecated in a future release of Biopython in favor of 'Bio.PDB.Polypeptide.protein_letters_3to1'.
  warnings.warn(
D:\Programming\anaconda3\envs\prot\Lib\site-packages\Bio\PDB\Polypeptide.py:144: BiopythonDeprecationWarning: 'three_to_one' will be deprecated in a future release of Biopython in favor of 'Bio.PDB.Polypeptide.protein_letters_3to1'.
  warnings.warn(
28/07/2024 18:41:36 logger: ============>chain 1 L and chain2 H
D:\Programming\anaconda3\envs\prot\Lib\site-packages\Bio\PDB\Polypeptide.py:144: BiopythonDeprecationWarning: 'three_to_one' will be deprecated in a future release of Biopython in favor of 'Bio.PDB.Polypeptide.protein_letters_3to1'.
  warnings.warn(
28/07/2024 18:41:46 logger: ============>chain 1 M and chain2 H
       res_fas_id1 res1  res_pdb_id1  ...       dist chain1  chain2
0                1    A            1  ...  33.363850      L       M
1                1    A            1  ...  27.019474      L       M
2                1    A            1  ...  29.652020      L       M
3                1    A            1  ...  28.374273      L       M
4                1    A            1  ...  29.460558      L       M
...            ...  ...          ...  ...        ...    ...     ...
74041          301    H          301  ...  61.372215      M       H
74042          301    H          301  ...  64.072952      M       H
74043          301    H          301  ...  62.540863      M       H
74044          301    H          301  ...  60.712658      M       H
74045          301    H          301  ...  65.942444      M       H

[227753 rows x 9 columns]
```

!!! tip

    All chains of a protein can be seen by `pp.str.chains`.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

chains = pp.str.chains(
    pdb_fp=to('data/pdb/complex/pdbtm/'),
    pdb_name='1aij',
)
print(chains)
```

:material-note-multiple-outline: Output
``` shell
['L', 'M', 'H']
```


## 2. Check interaction

It checks if a residue in the chain 1 has a minimum distance against all of residues in the chain 2. It stops the calculations of the minimum distance of each residue to each residue in the chain 2 when it detects a minimum distance of less than `thres`, 6 by default. `pp.dist.check_chain_complex` checks every two chains within a complex.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

pp.dist.check_chain_complex(
    pdb_fp=to('data/pdb/complex/pdbtm/'),
    prot_name='1aij',
    sv_fp=to('data/pdb/complex/pdbtm/'),
    thres=5.5,
)
```

:material-note-multiple-outline: Output
``` shell
28/07/2024 19:57:40 logger: =========>Protein PDB code: 1aij
D:\Programming\anaconda3\envs\prot\Lib\site-packages\Bio\PDB\Atom.py:232: PDBConstructionWarning: Could not assign element 'M' for Atom (name=MG) with given element 'M'
  warnings.warn(msg, PDBConstructionWarning)
28/07/2024 19:57:40 logger: =========>Protein chain 1: L
28/07/2024 19:57:40 logger: ============>Protein chain 2: M
28/07/2024 19:57:40 logger: ==================>residue 1 ID: 0
28/07/2024 19:57:40 logger: ==================>residue 0 and residue 252 in interaction
28/07/2024 19:57:40 logger: ============>Protein chain 2: H
28/07/2024 19:57:40 logger: ==================>residue 1 ID: 0
28/07/2024 19:57:40 logger: ==================>residue 0 and residue 31 in interaction
28/07/2024 19:57:40 logger: =========>Protein chain 1: M
28/07/2024 19:57:40 logger: ============>Protein chain 2: H
28/07/2024 19:57:40 logger: ==================>residue 1 ID: 0
28/07/2024 19:57:40 logger: ==================>residue 0 and residue 186 in interaction
28/07/2024 19:57:40 logger: =========>Protein chain 1: H
```

The results are saved to the `1aij.ccheck` file. It has the following content.
``` text
1aij	L
1aij	M
1aij	L
1aij	H
1aij	M
1aij	H
```

## 3. `dist` file

We can calculate the distance of a residue from a target protein away from residues from the rest of chains in the protein complex.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

pp.dist.complex_calc_all(
    pdb_fp=to('data/pdb/complex/pdbtm/'),
    prot_name='1aij',
    prot_chain='L',
    method='heavy',
    sv_fp=to('data/pdb/complex/pdbtm/'),
)
```

It returns a file `1aijL.dist`, having the following content.

:material-note-multiple-outline: Output
``` shell
       0  1    2         M          H
0      1  A    1  3.359507   3.820152
1      2  L    2  4.832898   2.840868
2      3  L    3  3.621853   3.450767
3      4  S    4  6.046180   2.841488
4      5  F    5  3.438960   3.688982
..   ... ..  ...       ...        ...
276  277  G  277  2.910155  43.177975
277  278  G  278  3.036153  40.448544
278  279  I  279  2.979592  36.381111
279  280  N  280  2.980944  38.922234
280  281  G  281  3.484761  42.785561

[281 rows x 5 columns]
```

The above columns correspond to different items in the following table.

| Column | Description      |
|:-------|:-----------------|
| 0      | fasta id 1       |
| 1      | residue 1        |
| 2      | pdb id 1         |
| M      | chain 1 distance |
| H      | chain 2 distance |


## 4. `pdbinter` file

We can calculate the distance of a residue from a target protein away from residues from the rest of chains in the protein complex.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

pp.dist.complex_calc_all(
    pdb_fp=to('data/pdb/complex/pdbtm/'),
    prot_name='1aij',
    prot_chain='L',
    method='heavy',
    sv_fp=to('data/pdb/complex/pdbtm/'),
)
```

It returns a file `1aijL.pdbinter`, having the following content.

:material-note-multiple-outline: Output
``` shell
L	1	A	M	253	R	3.3595068
L	2	L	M	253	R	4.832898
L	3	L	M	246	E	5.9632444
L	3	L	M	249	A	4.6895933
L	3	L	M	250	L	3.7801166
...
L	227	L	H	173	E	3.6427267
L	227	L	H	175	M	3.9549832
L	227	L	H	177	R	5.680695
L	228	G	H	173	E	5.4361606

[643 rows x 7 columns]
```

The above columns correspond to different items in the following table.

| Column  | Description |
|:--------|:------------|
| chain1  | chain 1     |
| pos1    | pdb id 1    |
| aa1     | residue 1   |
| chain2  | chain 2     |
| pos2    | pdb id 2    |
| aa2     | residue 2   |
| dist    | distance    |


## 5. Submission to the server

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

from pypropel.util.Reader import Reader as pfreader

df = pfreader().generic(df_fpn=to('data/ex/final.txt'), header=0)
prots = df.prot.unique()[2000:]

param_config = {
    'pdb_fp': '-fp',
    'pdb_fn': '-fn',
    'sv_fp': '-op',
}
value_config = {
    'tool_fp': '/path/to/python',
    'script_fpn': './Complex.py',
    'pdb_fp': '/path/to/protein complex files/',
    'sv_fp': '/path/to/save/results/',
}

for key, prot in enumerate(prots):
    order_list = [
        value_config['tool_fp'],
        value_config['script_fpn'],

        param_config['pdb_fp'], value_config['pdb_fp'],
        param_config['pdb_fn'], prot,
        param_config['sv_fp'], value_config['sv_fp'],
    ]
    pp.dist.cloud_check(
        order_list=order_list,
        job_fp='/path/to/save/job files/',
        job_fn=str(key),
        cpu=2,
        memory=10,
        method='script',
        submission_method='sbatch',
    )
```

!!! tip

    the `final.txt` file is obtained from running TMKit's `tmkit.collate` for a set of proteins.


``` shell
Unnamed: 0	prot	chain	prot_mark	rez	met	bio_name	head	desc	mthm	seq_aa	seq_len	nchain	met1
0	2a06	C	2a06C	2.1	x-ray diffraction	bovine cytochrome bc1 complex with stigmatellin bound	oxidoreductase	OXIDOREDUCTASE	8	NNAFIDLPAPSNISSWWNFGSLLGICLILQILTGLFLAMHYTSDTTTAFSSVTHICRDVNYGWIIRYMHANGASMFFICLYMHVGRGLYYGSYTFLETWNIGVILLLTVMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGTNLVEWIWGGFSVDKATLTRFFAFHFILPFIIMAIAMVHLLFLHETGSNNPTGISSDVDKIPFHPYYTIKDILGALLLILALMLLVLFAPDLLGDPDNYTPANPLNTPPHIKPEWYFLFAYAILRSIPNKLGGVLALAFSILILALIPLLHTSKQRSMMFRPLSQCLFWALVADLLTLTWIGGQPVEHPYITIGQLASVLYFLLILVLMPTAGTIENKLLKW	365.0	20.0	x-ray diffraction
4	2a06	P	2a06P	2.1		bovine cytochrome bc1 complex with stigmatellin bound	oxidoreductase	OXIDOREDUCTASE	8	LMKIVNNAFIDLPAPSNISSWWNFGSLLGICLILQILTGLFLAMHYTSDTTTAFSSVTHICRDVNYGWIIRYMHANGASMFFICLYMHVGRGLYYGSYTFLETWNIGVILLLTVMATAFMGYVLPWGQMSFWGATVITNLLSAIPYIGTNLVEWIWGGFSVDKATLTRFFAFHFILPFIIMAIAMVHLLFLHETGSNNPTGISSDVDKIPFHPYYTIKDILGALLLILALMLLVLFAPDLLGDPDNYTPANPLNTPPHIKPEWYFLFAYAILRSIPNKLGGVLALAFSILILALIPLLHTSKQRSMMFRPLSQCLFWALVADLLTLTWIGGQPVEHPYITIGQLASVLYFLLILVLMPTAGTIENKLLKW	370.0	20.0	x-ray diffraction
77	4a01	A	4a01A	2.35	x-ray diffraction	crystal structure of the h-translocating pyrophosphatase	hydrolase	HYDROLASE	16	GAAILPDLGTEILIPVCAVIGIAFALFQWLLVSKVKLSAVDHNVVVKCAEIQNAISEGATSFLFTEYKYVGIFMVAFAILIFLFLGSVEGFSTSPQACSYDKTKTCKPALATAIFSTVSFLLGGVTSLVSGFLGMKIATYANARTTLEARKGVGKAFITAFRSGAVMGFLLAANGLLVLYIAINLFKIYYGDDWGGLFEAITGYGLGGSSMALFGRVGGGIYTKAADVGADLVGKVERNIPEDDPRNPAVIADNVGDNVGDIAGMGSDLFGSYAESSCAALVVASISSFGLNHELTAMLYPLIVSSVGILVCLLTTLFATDFFEIKAVKEIEPALKKQLVISTVLMTIGVAVVSFVALPTSFTIFNFGVQKDVKSWQLFLCVAVGLWAGLIIGFVTEYYTSNAYSPVQDVADSCRTGAATNVIFGLALGYKSVIIPIFAIAISIFVSFTFAAMYGIAVAALGMLSTIATGLAIDAYGPISDNAGGIAEMAGMSHRIRERTDALDAAGNTTAAIGKGFAIGSAALVSLALFGAFVSRASITTVDVLTPKVFIGLIVGAMLPYWFSAMTMKSVGSAALKMVEEVRRQFNTIPGLMEGTAKPDYATCVKISTDASIKEMIPPGALVMLTPLVVGILFGVETLSGVLAGSLVSGVQIAISASNTGGAWDNAKKYIEAGASEHARSLGPKGSDCHKAAVIGDTIGDPLKDTSGPSLNILIKLMAVESLVFAPFFATHGGLLFKIF	740.0	2.0	x-ray diffraction
78	4a01	B	4a01B	2.35		crystal structure of the h-translocating pyrophosphatase	hydrolase	HYDROLASE	16	GAAILPDLGTEILIPVCAVIGIAFALFQWLLVSKVKLSAVDHNVVVKCAEIQNAISEGATSFLFTEYKYVGIFMVAFAILIFLFLGSVEGFSTSPQACSYDKTKTCKPALATAIFSTVSFLLGGVTSLVSGFLGMKIATYANARTTLEARKGVGKAFITAFRSGAVMGFLLAANGLLVLYIAINLFKIYYGDDWGGLFEAITGYGLGGSSMALFGRVGGGIYTKAADVGADLVGKVERNIPEDDPRNPAVIADNVGDNVGDIAGMGSDLFGSYAESSCAALVVASISSFGLNHELTAMLYPLIVSSVGILVCLLTTLFATDFFEIKAVKEIEPALKKQLVISTVLMTIGVAVVSFVALPTSFTIFNFGVQKDVKSWQLFLCVAVGLWAGLIIGFVTEYYTSNAYSPVQDVADSCRTGAATNVIFGLALGYKSVIIPIFAIAISIFVSFTFAAMYGIAVAALGMLSTIATGLAIDAYGPISDNAGGIAEMAGMSHRIRERTDALDAAGNTTAAIGKGFAIGSAALVSLALFGAFVSRASITTVDVLTPKVFIGLIVGAMLPYWFSAMTMKSVGSAALKMVEEVRRQFNTIPGLMEGTAKPDYATCVKISTDASIKEMIPPGALVMLTPLVVGILFGVETLSGVLAGSLVSGVQIAISASNTGGAWDNAKKYIEAGASEHARSLGPKGSDCHKAAVIGDTIGDPLKDTSGPSLNILIKLMAVESLVFAPFFATHGGLLFKIF	740.0	2.0	x-ray diffraction
79	7a0w	A	7a0wA	2.04	x-ray diffraction	structure of dimeric sodium proton antiporter nhaa, at ph 8.5, crystallized with chimeric fab antibodies	membrane protein	MEMBRANE PROTEIN	12	ASGGIILIIAAALAMLMANMGATSGWYHDFLETPVQLRVGALEINKNMLLWINDALMAVFFLLIGLEVKRELMQGSLASLRQAAFPVIAAIGGMIVPALLYLAFNYSDPVTREGWAIPAATDIAFALGVLALLGSRVPLALKIFLMALAIIDDLGAIVIIALFYTSDLSIVSLGVAAFAIAVLALLNLCGVRRTGVYILVGAVLWTAVLKSGVHATLAGVIVGFFIPLKEKHGRSPAKRLEHVLHPWVAYLILPLFAFANAGVSLQGVTIDGLTSMLPLGIIAGLLIGKPLGISLFCWLALRFKLAHLPQGTTYQQIMAVGILCGIGFTMSIFIASLAFGNVDPELINWAKLGILIGSLLSAVVGYSW	368.0	2.0	x-ray diffraction
```