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
{'AA': 0.0, 'AC': 0.0, 'AD': 0.0, 'AE': 0.0, 'AF': 0.0, 'AG': 0.0, 'AH': 0.0, 'AI': 0.0, 'AK': 0.0, 'AL': 0.0, 'AM': 0.0, 'AN': 0.0, 'AP': 0.0, 'AQ': 0.0, 'AR': 0.0, 'AS': 0.0, 'AT': 0.0, 'AV': 0.0, 'AW': 0.0, 'AY': 0.0, 'A-': 0.0, 'CA': 0.0, 'CC': 0.0, 'CD': 0.0, 'CE': 0.0, 'CF': 0.0, 'CG': 0.0, 'CH': 0.0, 'CI': 0.0, 'CK': 0.0, 'CL': 0.0, 'CM': 0.0, 'CN': 0.0, 'CP': 0.0, 'CQ': 0.0, 'CR': 0.0, 'CS': 0.0, 'CT': 0.0, 'CV': 0.0, 'CW': 0.0, 'CY': 0.0, 'C-': 0.0, 'DA': 0.0, 'DC': 0.0, 'DD': 0.0, 'DE': 0.0, 'DF': 0.0, 'DG': 0.0, 'DH': 0.0, 'DI': 0.0, 'DK': 0.0, 'DL': 0.0, 'DM': 0.0, 'DN': 0.0, 'DP': 0.0, 'DQ': 0.0, 'DR': 0.0, 'DS': 0.0, 'DT': 0.0, 'DV': 0.0, 'DW': 0.0, 'DY': 0.0, 'D-': 0.0, 'EA': 0.0, 'EC': 0.0, 'ED': 0.0, 'EE': 0.0, 'EF': 0.0, 'EG': 0.0, 'EH': 0.0, 'EI': 0.0, 'EK': 0.0, 'EL': 0.0, 'EM': 0.0, 'EN': 0.0, 'EP': 0.0, 'EQ': 0.0, 'ER': 0.0, 'ES': 0.0, 'ET': 0.0, 'EV': 0.0, 'EW': 0.0, 'EY': 0.0, 'E-': 0.0, 'FA': 0.0, 'FC': 0.0, 'FD': 0.0, 'FE': 0.0, 'FF': 0.0, 'FG': 0.0, 'FH': 0.0, 'FI': 0.0, 'FK': 0.0, 'FL': 0.0, 'FM': 0.0, 'FN': 0.0, 'FP': 0.0, 'FQ': 0.0, 'FR': 0.0, 'FS': 0.0, 'FT': 0.0, 'FV': 0.0, 'FW': 0.0, 'FY': 0.0, 'F-': 0.0, 'GA': 0.0, 'GC': 0.0, 'GD': 0.0, 'GE': 0.0, 'GF': 0.0, 'GG': 0.0, 'GH': 0.0, 'GI': 0.0, 'GK': 0.0, 'GL': 0.0, 'GM': 0.0, 'GN': 0.0, 'GP': 0.0, 'GQ': 0.0, 'GR': 0.0, 'GS': 0.0, 'GT': 0.0, 'GV': 0.0, 'GW': 0.0, 'GY': 0.0, 'G-': 0.0, 'HA': 0.0, 'HC': 0.0, 'HD': 0.0, 'HE': 0.0, 'HF': 0.0, 'HG': 0.0, 'HH': 0.0, 'HI': 0.0, 'HK': 0.0, 'HL': 0.0, 'HM': 0.0, 'HN': 0.0, 'HP': 0.0, 'HQ': 0.0, 'HR': 0.0, 'HS': 0.0, 'HT': 0.0, 'HV': 0.0, 'HW': 0.0, 'HY': 0.0, 'H-': 0.0, 'IA': 0.0, 'IC': 0.0, 'ID': 0.0, 'IE': 0.0, 'IF': 0.0, 'IG': 0.0, 'IH': 0.0, 'II': 0.0, 'IK': 0.0, 'IL': 0.0002472035102898461, 'IM': 0.0, 'IN': 0.0, 'IP': 0.0, 'IQ': 0.0, 'IR': 0.0, 'IS': 0.0, 'IT': 0.0, 'IV': 0.0, 'IW': 0.0, 'IY': 0.0, 'I-': 0.0, 'KA': 0.0, 'KC': 0.0, 'KD': 0.0, 'KE': 0.0, 'KF': 0.0, 'KG': 0.0, 'KH': 0.0, 'KI': 0.0, 'KK': 0.0, 'KL': 0.0, 'KM': 0.0, 'KN': 0.0, 'KP': 0.0, 'KQ': 0.0, 'KR': 0.0, 'KS': 0.0, 'KT': 0.0, 'KV': 0.0, 'KW': 0.0, 'KY': 0.0, 'K-': 0.0, 'LA': 0.0, 'LC': 0.0, 'LD': 0.0, 'LE': 0.0, 'LF': 0.0, 'LG': 0.0, 'LH': 0.0, 'LI': 0.0, 'LK': 0.0, 'LL': 0.025585563314999074, 'LM': 0.0, 'LN': 0.0, 'LP': 0.0, 'LQ': 0.0, 'LR': 0.0, 'LS': 0.0, 'LT': 0.0, 'LV': 0.0, 'LW': 0.0, 'LY': 0.0, 'L-': 0.0, 'MA': 0.0, 'MC': 0.0, 'MD': 0.0, 'ME': 0.0, 'MF': 0.0, 'MG': 0.0, 'MH': 0.0, 'MI': 0.0, 'MK': 0.0, 'ML': 0.023916939620542612, 'MM': 0.0008652122860144614, 'MN': 0.0, 'MP': 0.0, 'MQ': 0.0, 'MR': 0.0, 'MS': 0.0, 'MT': 0.0, 'MV': 0.0, 'MW': 0.0, 'MY': 0.0, 'M-': 0.0, 'NA': 0.0, 'NC': 0.0, 'ND': 0.0, 'NE': 0.0, 'NF': 0.0, 'NG': 0.0, 'NH': 0.0, 'NI': 0.0, 'NK': 0.0, 'NL': 0.0, 'NM': 0.0, 'NN': 0.0, 'NP': 0.0, 'NQ': 0.0, 'NR': 0.0, 'NS': 0.0, 'NT': 0.0, 'NV': 0.0, 'NW': 0.0, 'NY': 0.0, 'N-': 0.0, 'PA': 0.0, 'PC': 0.0, 'PD': 0.0, 'PE': 0.0, 'PF': 0.0, 'PG': 0.0, 'PH': 0.0, 'PI': 0.0, 'PK': 0.0, 'PL': 0.0, 'PM': 0.0, 'PN': 0.0, 'PP': 0.0, 'PQ': 0.0, 'PR': 0.0, 'PS': 0.0, 'PT': 0.0, 'PV': 0.0, 'PW': 0.0, 'PY': 0.0, 'P-': 0.0, 'QA': 0.0, 'QC': 0.0, 'QD': 0.0, 'QE': 0.0, 'QF': 0.0, 'QG': 0.0, 'QH': 0.0, 'QI': 0.0, 'QK': 0.0, 'QL': 0.002101229837463692, 'QM': 0.0, 'QN': 0.0, 'QP': 0.0, 'QQ': 0.0, 'QR': 0.0, 'QS': 0.0, 'QT': 0.0, 'QV': 0.0, 'QW': 0.0, 'QY': 0.0, 'Q-': 0.0, 'RA': 0.0, 'RC': 0.0, 'RD': 0.0, 'RE': 0.0, 'RF': 0.0, 'RG': 0.0, 'RH': 0.0, 'RI': 0.0, 'RK': 0.0, 'RL': 0.0, 'RM': 0.0, 'RN': 0.0, 'RP': 0.0, 'RQ': 0.0, 'RR': 0.0, 'RS': 0.0, 'RT': 0.0, 'RV': 0.0, 'RW': 0.0, 'RY': 0.0, 'R-': 0.0, 'SA': 0.0, 'SC': 0.0, 'SD': 0.0, 'SE': 0.0, 'SF': 0.0, 'SG': 0.0, 'SH': 0.0, 'SI': 0.0, 'SK': 0.0, 'SL': 0.0, 'SM': 0.0, 'SN': 0.0, 'SP': 0.0, 'SQ': 0.0, 'SR': 0.0, 'SS': 0.0, 'ST': 0.0, 'SV': 0.0, 'SW': 0.0, 'SY': 0.0, 'S-': 0.0, 'TA': 0.0, 'TC': 0.0, 'TD': 0.0, 'TE': 0.0, 'TF': 0.0, 'TG': 0.0, 'TH': 0.0, 'TI': 0.0, 'TK': 0.0, 'TL': 0.00018540263271738457, 'TM': 0.0, 'TN': 0.0, 'TP': 0.0, 'TQ': 0.0, 'TR': 0.0, 'TS': 0.0, 'TT': 0.0, 'TV': 0.0, 'TW': 0.0, 'TY': 0.0, 'T-': 0.0, 'VA': 0.0, 'VC': 0.0, 'VD': 0.0, 'VE': 0.0, 'VF': 0.0, 'VG': 0.0, 'VH': 0.0, 'VI': 0.0, 'VK': 0.0, 'VL': 0.0005562078981521537, 'VM': 0.0, 'VN': 0.0, 'VP': 0.0, 'VQ': 0.0, 'VR': 0.0, 'VS': 0.0, 'VT': 0.0, 'VV': 0.0, 'VW': 0.0, 'VY': 0.0, 'V-': 0.0, 'WA': 0.0, 'WC': 0.0, 'WD': 0.0, 'WE': 0.0, 'WF': 0.0, 'WG': 0.0, 'WH': 0.0, 'WI': 0.0, 'WK': 0.0, 'WL': 0.0, 'WM': 0.0, 'WN': 0.0, 'WP': 0.0, 'WQ': 0.0, 'WR': 0.0, 'WS': 0.0, 'WT': 0.0, 'WV': 0.0, 'WW': 0.0, 'WY': 0.0, 'W-': 0.0, 'YA': 0.0, 'YC': 0.0, 'YD': 0.0, 'YE': 0.0, 'YF': 0.0, 'YG': 0.0, 'YH': 0.0, 'YI': 0.0, 'YK': 0.0, 'YL': 0.0, 'YM': 0.0, 'YN': 0.0, 'YP': 0.0, 'YQ': 0.0, 'YR': 0.0, 'YS': 0.0, 'YT': 0.0, 'YV': 0.0, 'YW': 0.0, 'YY': 0.0, 'Y-': 0.0, '-A': 0.0, '-C': 0.0, '-D': 0.0, '-E': 0.0, '-F': 0.0, '-G': 0.0, '-H': 0.0, '-I': 0.0, '-K': 0.0, '-L': 0.0, '-M': 0.0, '-N': 0.0, '-P': 0.0, '-Q': 0.0, '-R': 0.0, '-S': 0.0, '-T': 0.0, '-V': 0.0, '-W': 0.0, '-Y': 0.0, '--': 0.9465422408998208}
```