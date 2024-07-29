We use 23 indicators to describe physicochemical properties of amino acids as shown in Table 1. You can use `pp.fpsite.property` to get access to them.

:material-language-python: Python
``` py linenums="1"
import pypropel as pp

property = 'positive'

scale = pp.fpsite.property(property)

print(scale)
```

:material-note-multiple-outline: Output
``` shell
{'A': 1.0, 'C': 1.0, 'D': 1.0, 'E': 1.0, 'F': 1.0, 'G': 1.0, 'H': 0.5, 'I': 1.0, 'K': 0.0, 'L': 1.0, 'M': 1.0, 'N': 1.0, 'P': 1.0, 'Q': 1.0, 'R': 0.0, 'S': 1.0, 'T': 1.0, 'V': 1.0, 'W': 1.0, 'Y': 1.0}
```

In fact, `property` can be any of 23 amino acid properties, including

**Table 1**. amino acid properties

| Property        | Kind         | Citation                  | Value     |
|:----------------|:-------------|:--------------------------|-----------|
| positive        | `Russell`    | Betts & Russell[^1]       | discrete  |
| negative        | `Russell`    | Betts & Russell[^1]       | discrete  |
| charged         | `Russell`    | Betts & Russell[^1]       | discrete  |
| polar           | `Russell`    | Betts & Russell[^1]       | discrete  |
| aliphatic       | `Russell`    | Betts & Russell[^1]       | discrete  |
| aromatic        | `Russell`    | Betts & Russell[^1]       | discrete  |
| hydrophobic     | `Russell`    | Betts & Russell[^1]       | discrete  |
| small           | `Russell`    | Betts & Russell[^1]       | discrete  |
| active          | `Russell`    | Betts & Russell[^1]       | continuous |
| weight          | `Taylor`     | Lundblad & Macdonald[^2]  | continuous |
| pI              | `Taylor`     | Lundblad & Macdonald[^2]  | continuous |
| solubility      | `Taylor`     | Lundblad & Macdonald[^2]  | continuous |
| tm              | `Taylor`     | Lundblad & Macdonald[^2]  | continuous |
| pka             | `Taylor`     | Lundblad & Macdonald[^2]  | continuous |
| pkb             | `Taylor`     | Lundblad & Macdonald[^2]  | continuous |
| hydrophilicity  | `Hopp`       | [^3]                      | discrete  |
| hydrophobicity  | `Argos`      | Argos et al.[^4]          | continuous |
| fet             | `Argos`      | Argos et al.[^4]          | continuous |
| hydration       | `Argos`      | Argos et al.[^4]          | continuous |
| signal          | `Argos`      | Argos et al.[^4]          | continuous |
| volume          | `Grantham`   | Grantham [^5]             | continuous |
| polarity        | `Grantham`   | Grantham [^5]             | continuous |
| composition     | `Grantham`   | Grantham [^5]             | continuous |

!!! note

    pI: pH at the isoelectric point.
    
    Solubility in water in units of grams of compound per kilogram of water.
    
    tm: Melting point.
    
    pKa: Negative of the logarithm of the acid dissociation constants for the COOH and NH2 groups (and, in some cases, other groups) in the molecule (at 25°C)

[^1]: Betts, M.J. and Russell, R.B. (2003). Amino Acid Properties and Consequences of Substitutions. In Bioinformatics for Geneticists (eds M.R. Barnes and I.C. Gray). https://doi.org/10.1002/0470867302.ch14

[^2]: Lundblad, R.L., & Macdonald, F. (Eds.). (2018). Handbook of Biochemistry and Molecular Biology (5th ed.). CRC Press. https://doi.org/10.1201/b21846

[^3]: Hopp TP, Woods KR. Prediction of protein antigenic determinants from amino acid sequences. Proc Natl Acad Sci U S A. 1981 Jun;78(6):3824-8. doi: 10.1073/pnas.78.6.3824.

[^4]: ARGOS, P., RAO, J.K.M. and HARGRAVE, P.A. (1982), Structural Prediction of Membrane-Bound Proteins. European Journal of Biochemistry, 128: 565-575. https://doi.org/10.1111/j.1432-1033.1982.tb07002.x

[^5]: R. Grantham ,Amino Acid Difference Formula to Help Explain Protein Evolution.Science185,862-864(1974).DOI:10.1126/science.185.4154.862


