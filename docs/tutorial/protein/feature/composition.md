## Protein-level
In PyPropel, the function `fpseq` governs the extraction of features from the whole proteins. At the protein level, PyPropel considers compositions of amino acids across the whole protein sequence. It can be done below.

### AAC
AAC stands for the amino acid composition.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

aac = pp.fpseq.composition(
    seq=seq,
    mode='aac',
)
```

:material-note-multiple-outline: Output
``` shell
{'A': 0.119403, 'C': 0.074627, 'D': 0.149254, 'E': 0.089552, 'F': 0.059701, 'G': 0.089552, 'H': 0.0, 'I': 0.0, 'K': 0.014925, 'L': 0.029851, 'M': 0.044776, 'N': 0.029851, 'P': 0.029851, 'Q': 0.029851, 'R': 0.014925, 'S': 0.119403, 'T': 0.014925, 'V': 0.044776, 'W': 0.029851, 'Y': 0.014925}
```

### DAC
DAC stands for the di-amino acid composition.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

dac = pp.msa.composition(
    seq=seq,
    mode='dac',
)
```

:material-note-multiple-outline: Output
``` shell
[['AA', 0.015152], ['AC', 0.0], ['AD', 0.075758], ..., ['YW', 0.0], ['YY', 0.0]]
```

### TAC
TAC stands for the tri-amino acid composition.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

tac = pp.msa.composition(
    seq=seq,
    mode='tac',
)
```

:material-note-multiple-outline: Output
``` shell
[['AAA', 0.0], ['AAC', 0.0], ['AAD', 0.015385], ...]
```

### QAC
QAC stands for the qua-amino acid composition.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

qac = pp.msa.composition(
    seq=seq,
    mode='qac',
)
```

:material-note-multiple-outline: Output
``` shell
[['AAAA', 0.0], ['AAAC', 0.0], ['AAAD', 0.0], ...]
```

### CKSNAP
CKSNAP measures the average frequencies of every two amino acids that are `k` spaced.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

cksnap = pp.msa.composition(
    seq=seq,
    k_spaced=3,
    mode='cksnap',
)
```

:material-note-multiple-outline: Output
``` shell
{'AA': 0.047619, 'AC': 0, 'AD': 0, 'AE': 0.015873, 'AF': 0, 'AG': 0.015873, 'AH': 0, 'AI': 0, 'AK': 0, 'AL': 0, 'AM': 0, 'AN': 0, 'AP': 0, 'AQ': 0, 'AR': 0, 'AS': 0.031746, 'AT': 0, 'AV': 0, 'AW': 0, 'AY': 0, 'CA': 0, 'CC': 0.015873, 'CD': 0.015873, 'CE': 0.015873, 'CF': 0, 'CG': 0, 'CH': 0, 'CI': 0, 'CK': 0.015873, 'CL': 0, 'CM': 0, 'CN': 0, 'CP': 0, 'CQ': 0, 'CR': 0, 'CS': 0, 'CT': 0, 'CV': 0.015873, 'CW': 0, 'CY': 0, 'DA': 0.015873, 'DC': 0, 'DD': 0.031746, 'DE': 0.015873, 'DF': 0, 'DG': 0, 'DH': 0, 'DI': 0, 'DK': 0, 'DL': 0.015873, 'DM': 0, 'DN': 0, 'DP': 0, 'DQ': 0, 'DR': 0, 'DS': 0.031746, 'DT': 0, 'DV': 0.015873, 'DW': 0.015873, 'DY': 0, 'EA': 0, 'EC': 0, 'ED': 0.015873, 'EE': 0, 'EF': 0.031746, 'EG': 0, 'EH': 0, 'EI': 0, 'EK': 0, 'EL': 0.015873, 'EM': 0, 'EN': 0, 'EP': 0, 'EQ': 0.015873, 'ER': 0, 'ES': 0, 'ET': 0, 'EV': 0, 'EW': 0, 'EY': 0, 'FA': 0.015873, 'FC': 0, 'FD': 0.031746, 'FE': 0, 'FF': 0, 'FG': 0, 'FH': 0, 'FI': 0, 'FK': 0, 'FL': 0, 'FM': 0, 'FN': 0.015873, 'FP': 0, 'FQ': 0, 'FR': 0, 'FS': 0, 'FT': 0, 'FV': 0, 'FW': 0, 'FY': 0, 'GA': 0, 'GC': 0, 'GD': 0, 'GE': 0, 'GF': 0, 'GG': 0.063492, 'GH': 0, 'GI': 0, 'GK': 0, 'GL': 0, 'GM': 0.015873, 'GN': 0.015873, 'GP': 0, 'GQ': 0, 'GR': 0, 'GS': 0, 'GT': 0, 'GV': 0, 'GW': 0, 'GY': 0, 'HA': 0, 'HC': 0, 'HD': 0, 'HE': 0, 'HF': 0, 'HG': 0, 'HH': 0, 'HI': 0, 'HK': 0, 'HL': 0, 'HM': 0, 'HN': 0, 'HP': 0, 'HQ': 0, 'HR': 0, 'HS': 0, 'HT': 0, 'HV': 0, 'HW': 0, 'HY': 0, 'IA': 0, 'IC': 0, 'ID': 0, 'IE': 0, 'IF': 0, 'IG': 0, 'IH': 0, 'II': 0, 'IK': 0, 'IL': 0, 'IM': 0, 'IN': 0, 'IP': 0, 'IQ': 0, 'IR': 0, 'IS': 0, 'IT': 0, 'IV': 0, 'IW': 0, 'IY': 0, 'KA': 0.015873, 'KC': 0, 'KD': 0, 'KE': 0, 'KF': 0, 'KG': 0, 'KH': 0, 'KI': 0, 'KK': 0, 'KL': 0, 'KM': 0, 'KN': 0, 'KP': 0, 'KQ': 0, 'KR': 0, 'KS': 0, 'KT': 0, 'KV': 0, 'KW': 0, 'KY': 0, 'LA': 0, 'LC': 0, 'LD': 0, 'LE': 0, 'LF': 0, 'LG': 0, 'LH': 0, 'LI': 0, 'LK': 0, 'LL': 0, 'LM': 0, 'LN': 0, 'LP': 0, 'LQ': 0, 'LR': 0, 'LS': 0.031746, 'LT': 0, 'LV': 0, 'LW': 0, 'LY': 0, 'MA': 0, 'MC': 0, 'MD': 0, 'ME': 0, 'MF': 0, 'MG': 0, 'MH': 0, 'MI': 0, 'MK': 0, 'ML': 0, 'MM': 0.015873, 'MN': 0, 'MP': 0, 'MQ': 0, 'MR': 0, 'MS': 0, 'MT': 0, 'MV': 0, 'MW': 0.015873, 'MY': 0.015873, 'NA': 0, 'NC': 0, 'ND': 0, 'NE': 0, 'NF': 0.015873, 'NG': 0, 'NH': 0, 'NI': 0, 'NK': 0, 'NL': 0, 'NM': 0.015873, 'NN': 0, 'NP': 0, 'NQ': 0, 'NR': 0, 'NS': 0, 'NT': 0, 'NV': 0, 'NW': 0, 'NY': 0, 'PA': 0, 'PC': 0.015873, 'PD': 0, 'PE': 0, 'PF': 0, 'PG': 0, 'PH': 0, 'PI': 0, 'PK': 0, 'PL': 0, 'PM': 0, 'PN': 0, 'PP': 0, 'PQ': 0.015873, 'PR': 0, 'PS': 0, 'PT': 0, 'PV': 0, 'PW': 0, 'PY': 0, 'QA': 0, 'QC': 0.031746, 'QD': 0, 'QE': 0, 'QF': 0, 'QG': 0, 'QH': 0, 'QI': 0, 'QK': 0, 'QL': 0, 'QM': 0, 'QN': 0, 'QP': 0, 'QQ': 0, 'QR': 0, 'QS': 0, 'QT': 0, 'QV': 0, 'QW': 0, 'QY': 0, 'RA': 0, 'RC': 0, 'RD': 0, 'RE': 0, 'RF': 0, 'RG': 0, 'RH': 0, 'RI': 0, 'RK': 0, 'RL': 0, 'RM': 0, 'RN': 0, 'RP': 0, 'RQ': 0, 'RR': 0, 'RS': 0.015873, 'RT': 0, 'RV': 0, 'RW': 0, 'RY': 0, 'SA': 0.015873, 'SC': 0, 'SD': 0.015873, 'SE': 0.031746, 'SF': 0, 'SG': 0, 'SH': 0, 'SI': 0, 'SK': 0, 'SL': 0, 'SM': 0, 'SN': 0, 'SP': 0.015873, 'SQ': 0, 'SR': 0, 'SS': 0.015873, 'ST': 0, 'SV': 0.015873, 'SW': 0, 'SY': 0, 'TA': 0, 'TC': 0, 'TD': 0, 'TE': 0, 'TF': 0, 'TG': 0, 'TH': 0, 'TI': 0, 'TK': 0, 'TL': 0, 'TM': 0, 'TN': 0, 'TP': 0.015873, 'TQ': 0, 'TR': 0, 'TS': 0, 'TT': 0, 'TV': 0, 'TW': 0, 'TY': 0, 'VA': 0, 'VC': 0, 'VD': 0, 'VE': 0.015873, 'VF': 0.015873, 'VG': 0, 'VH': 0, 'VI': 0, 'VK': 0, 'VL': 0, 'VM': 0, 'VN': 0, 'VP': 0, 'VQ': 0, 'VR': 0, 'VS': 0, 'VT': 0.015873, 'VV': 0, 'VW': 0, 'VY': 0, 'WA': 0, 'WC': 0, 'WD': 0.015873, 'WE': 0, 'WF': 0, 'WG': 0, 'WH': 0, 'WI': 0, 'WK': 0, 'WL': 0, 'WM': 0, 'WN': 0, 'WP': 0, 'WQ': 0, 'WR': 0.015873, 'WS': 0, 'WT': 0, 'WV': 0, 'WW': 0, 'WY': 0, 'YA': 0, 'YC': 0, 'YD': 0.015873, 'YE': 0, 'YF': 0, 'YG': 0, 'YH': 0, 'YI': 0, 'YK': 0, 'YL': 0, 'YM': 0, 'YN': 0, 'YP': 0, 'YQ': 0, 'YR': 0, 'YS': 0, 'YT': 0, 'YV': 0, 'YW': 0, 'YY': 0}
```

### AVEANF
The AVEANF feature reflects the information about positions and frequencies of amino acids in a protein sequence. For details, please refer to [^1]

[^1]: Zhen Chen, Pei Zhao, Fuyi Li, Tatiana T Marquez-Lago, André Leier, Jerico Revote, Yan Zhu, David R Powell, Tatsuya Akutsu, Geoffrey I Webb, Kuo-Chen Chou, A Ian Smith, Roger J Daly, Jian Li, Jiangning Song, iLearn: an integrated platform and meta-learner for feature engineering, machine-learning analysis and modeling of DNA, RNA and protein sequence data, Briefings in Bioinformatics, Volume 21, Issue 3, May 2020, Pages 1047–1057, https://doi.org/10.1093/bib/bbz041

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

aveanf = pp.msa.composition(
    seq=seq,
    mode='aveanf',
)
```

:material-note-multiple-outline: Output
``` shell
{'A': 0.22597633674152806, 'C': 0.1652824858757062, 'D': 0.17214478688616622, 'E': 0.09744508222164794, 'F': 0.0459034244741305, 'G': 0.4204055204055204, 'H': 0, 'I': 0, 'K': 0.045454545454545456, 'L': 0.03737373737373738, 'M': 0.10492898913951544, 'N': 0.04793028322440087, 'P': 0.055322128851540614, 'Q': 0.059848484848484845, 'R': 0.017857142857142856, 'S': 0.08590396237795747, 'T': 0.1, 'V': 0.10087719298245613, 'W': 0.04096989966555184, 'Y': 0.04}
```

## MSA-level

First, we can read a MSA of protein `1aijL` as follows.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

msa = pp.msa.read(
    msa_fpn=to('data/msa/aln/1aijL.aln'),
)
```

### AAC

Then, the AAC of MSA is calculated below.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

aac = pp.msa.composition(
    msa=msa,
    mode='aac',
)
```

:material-note-multiple-outline: Output
``` shell
{'A': 0.05826085292688736, 'C': 0.004700385606685579, 'D': 0.008600878716107662, 'E': 0.015672130729309737, 'F': 0.05294113895278523, 'G': 0.06534156201388167, 'H': 0.019759785927038455, 'I': 0.04494023459261235, 'K': 0.0027040633087310126, 'L': 0.05968535215833517, 'M': 0.022540825417799226, 'N': 0.01639108826946766, 'P': 0.03155605592517563, 'Q': 0.011760641022454832, 'R': 0.017077715813173088, 'S': 0.03935616241622517, 'T': 0.02184298134471232, 'V': 0.040838283818220966, 'W': 0.01871466930702302, 'Y': 0.0222107075628659, '-': 0.42443809036607894}
```

## Site-level
There is no method to calculate compositions of amino acids because no accurate einformation can be reflected from a site itself in a protein sequence given only.