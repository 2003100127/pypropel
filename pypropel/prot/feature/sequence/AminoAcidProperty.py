__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import time
from pypropel.util.normalization.Standardize import Standardize


class AminoAcidProperty:

    def __init__(self, ):
        self.sstandardize = Standardize()

    def positive(
            self,
            kind='Russell',
            standardize=True,
    ):
        """
        Positive

        References
        ----------
        Betts, M.J. and Russell, R.B. (2003). Amino Acid Properties and Consequences
        of Substitutions. In Bioinformatics for Geneticists
        (eds M.R. Barnes and I.C. Gray). https://doi.org/10.1002/0470867302.ch14

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Russell':
            weight_dict = {
                'A': 3, 'C': 3, 'D': 3, 'E': 3,
                'F': 3, 'G': 3, 'H': 2, 'I': 3,
                'K': 1, 'L': 3, 'M': 3, 'N': 3,
                'P': 3, 'Q': 3, 'R': 1, 'S': 3,
                'T': 3, 'V': 3, 'W': 3, 'Y': 3
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def negative(
            self,
            kind='Russell',
            standardize=True,
    ):
        """
        Negative

        References
        ----------
        Betts, M.J. and Russell, R.B. (2003). Amino Acid Properties and Consequences
        of Substitutions. In Bioinformatics for Geneticists
        (eds M.R. Barnes and I.C. Gray). https://doi.org/10.1002/0470867302.ch14

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Russell':
            weight_dict = {
                'A': 2, 'C': 2, 'D': 1, 'E': 1,
                'F': 2, 'G': 2, 'H': 2, 'I': 2,
                'K': 2, 'L': 2, 'M': 2, 'N': 2,
                'P': 2, 'Q': 2, 'R': 2, 'S': 2,
                'T': 2, 'V': 2, 'W': 2, 'Y': 2
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def charged(
            self,
            kind='Russell',
            standardize=True,
    ):
        """
        Charged

        References
        ----------
        Betts, M.J. and Russell, R.B. (2003). Amino Acid Properties and Consequences
        of Substitutions. In Bioinformatics for Geneticists
        (eds M.R. Barnes and I.C. Gray). https://doi.org/10.1002/0470867302.ch14

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Russell':
            weight_dict = {
                'A': 4, 'C': 4, 'D': 1, 'E': 1,
                'F': 4, 'G': 4, 'H': 3, 'I': 4,
                'K': 2, 'L': 4, 'M': 4, 'N': 4,
                'P': 4, 'Q': 4, 'R': 2, 'S': 4,
                'T': 4, 'V': 4, 'W': 4, 'Y': 4
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def polar(
            self,
            kind='Russell',
            standardize=True,
    ):
        """
        Polar

        References
        ----------
        Betts, M.J. and Russell, R.B. (2003). Amino Acid Properties and Consequences
        of Substitutions. In Bioinformatics for Geneticists
        (eds M.R. Barnes and I.C. Gray). https://doi.org/10.1002/0470867302.ch14

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Russell':
            weight_dict = {
                'A': 2, 'C': 3, 'D': 1, 'E': 1,
                'F': 3, 'G': 2, 'H': 2, 'I': 3,
                'K': 1, 'L': 3, 'M': 3, 'N': 1,
                'P': 2, 'Q': 1, 'R': 1, 'S': 2,
                'T': 2, 'V': 3, 'W': 3, 'Y': 2
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def aliphatic(
            self,
            kind='Russell',
            standardize=True,
    ):
        """
        Aliphatic

        References
        ----------
        Betts, M.J. and Russell, R.B. (2003). Amino Acid Properties and Consequences
        of Substitutions. In Bioinformatics for Geneticists
        (eds M.R. Barnes and I.C. Gray). https://doi.org/10.1002/0470867302.ch14

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Russell':
            weight_dict = {
                'A': 1, 'C': 2, 'D': 2, 'E': 2,
                'F': 2, 'G': 2, 'H': 2, 'I': 1,
                'K': 2, 'L': 1, 'M': 1, 'N': 2,
                'P': 1, 'Q': 2, 'R': 2, 'S': 2,
                'T': 2, 'V': 1, 'W': 2, 'Y': 2
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def aromatic(
            self,
            kind='Russell',
            standardize=True,
    ):
        """
        Aromatic

        References
        ----------
        Betts, M.J. and Russell, R.B. (2003). Amino Acid Properties and Consequences
        of Substitutions. In Bioinformatics for Geneticists
        (eds M.R. Barnes and I.C. Gray). https://doi.org/10.1002/0470867302.ch14

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Russell':
            weight_dict = {
                'A': 2, 'C': 2, 'D': 2, 'E': 2,
                'F': 1, 'G': 2, 'H': 1, 'I': 2,
                'K': 2, 'L': 2, 'M': 2, 'N': 2,
                'P': 2, 'Q': 2, 'R': 2, 'S': 2,
                'T': 2, 'V': 2, 'W': 1, 'Y': 1
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def hydrophobic(
            self,
            kind='Russell',
            standardize=True,
    ):
        """
        Hydrophobic

        References
        ----------
        Betts, M.J. and Russell, R.B. (2003). Amino Acid Properties and Consequences
        of Substitutions. In Bioinformatics for Geneticists
        (eds M.R. Barnes and I.C. Gray). https://doi.org/10.1002/0470867302.ch14

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Russell':
            weight_dict = {
                'A': 2, 'C': 1, 'D': 4, 'E': 4,
                'F': 1, 'G': 2, 'H': 2, 'I': 1,
                'K': 3, 'L': 1, 'M': 1, 'N': 4,
                'P': 2, 'Q': 4, 'R': 3, 'S': 2,
                'T': 2, 'V': 1, 'W': 1, 'Y': 2
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def small(
            self,
            kind='Russell',
            standardize=True,
    ):
        """
        Small

        References
        ----------
        Betts, M.J. and Russell, R.B. (2003). Amino Acid Properties and Consequences
        of Substitutions. In Bioinformatics for Geneticists
        (eds M.R. Barnes and I.C. Gray). https://doi.org/10.1002/0470867302.ch14

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Russell':
            weight_dict = {
                'A': 1, 'C': 3, 'D': 2, 'E': 3,
                'F': 3, 'G': 1, 'H': 3, 'I': 3,
                'K': 3, 'L': 3, 'M': 3, 'N': 2,
                'P': 1, 'Q': 3, 'R': 3, 'S': 1,
                'T': 2, 'V': 3, 'W': 3, 'Y': 3
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def active(
            self,
            kind='Russell',
            standardize=True,
    ):
        """
        Active

        References
        ----------
        Betts, M.J. and Russell, R.B. (2003). Amino Acid Properties and Consequences
        of Substitutions. In Bioinformatics for Geneticists
        (eds M.R. Barnes and I.C. Gray). https://doi.org/10.1002/0470867302.ch14

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Russell':
            weight_dict ={
                'A': 0.025, 'C': 0.21, 'D': 0.045, 'E': 0.05,
                'F': -0.12, 'G': -0.07, 'H': 0.36, 'I': -0.005,
                'K': 0.1, 'L': -0.18, 'M': 0.025, 'N': 0.08,
                'P': -0.2, 'Q': 0.05, 'R': 0.055, 'S': 0.13,
                'T': 0.1, 'V': -0.06, 'W': -0.14, 'Y': -0.04
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def weight(
            self,
            kind='Taylor',
            standardize=True,
    ):
        """
        Molecular weight

        References
        ----------
        Lundblad, R.L., & Macdonald, F. (Eds.). (2018). Handbook of Biochemistry
         and Molecular Biology (5th ed.). CRC Press. https://doi.org/10.1201/b21846

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Taylor':
            weight_dict = {
                'A': 89.09, 'C': 121.15, 'D': 133.10,
                'E': 147.13, 'F': 165.19, 'G': 75.07,
                'H': 155.16, 'I': 131.17, 'K': 146.19,
                'L': 131.17, 'M': 149.21, 'N': 132.12,
                'P': 115.13, 'Q': 146.15, 'R': 174.20,
                'S': 105.09, 'T': 119.12, 'V': 117.15,
                'W': 204.24, 'Y': 181.19
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def pI(
            self,
            kind='Taylor',
            standardize=True,
    ):
        """
        pI: pH at the isoelectric point

        References
        ----------
        Lundblad, R.L., & Macdonald, F. (Eds.). (2018). Handbook of Biochemistry
         and Molecular Biology (5th ed.). CRC Press. https://doi.org/10.1201/b21846

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Taylor':
            weight_dict = {
                'A': 6.00, 'C': 5.07, 'D': 2.77, 'E': 3.22,
                'F': 5.48, 'G': 5.97, 'H': 7.59, 'I': 6.02,
                'K': 9.74, 'L': 5.98, 'M': 5.74, 'N': 5.41,
                'P': 6.30, 'Q': 5.65, 'R': 10.76, 'S': 5.68,
                'T': 5.60, 'V': 5.96, 'W': 5.89, 'Y': 5.66
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def solubility(
            self,
            kind='Taylor',
            standardize=True,
    ):
        """
        Solubility in water in units of grams of
            compound per kilogram of water.

        References
        ----------
        Lundblad, R.L., & Macdonald, F. (Eds.). (2018). Handbook of Biochemistry
         and Molecular Biology (5th ed.). CRC Press. https://doi.org/10.1201/b21846

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Taylor':
            weight_dict = {
                'A': 166.9, 'C': 1623, 'D': 5.04, 'E': 8.6,
                'F': 27.9, 'G': 250.2, 'H': 43.5, 'I': 34.2,
                'K': 5.8, 'L': 22.0, 'M': 56, 'N': 25.1,
                'P': 1622, 'Q': 42.5, 'R': 182.6, 'S': 250,
                'T': 98.1, 'V': 88, 'W': 13.2, 'Y': 0.46
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def tm(
            self,
            kind='Taylor',
            standardize=True,
    ):
        """
        tm: Melting point

        References
        ----------
        Lundblad, R.L., & Macdonald, F. (Eds.). (2018). Handbook of Biochemistry
         and Molecular Biology (5th ed.). CRC Press. https://doi.org/10.1201/b21846

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Taylor':
            weight_dict = {
                'A': 297, 'C': 240, 'D': 270, 'E': 160,
                'F': 283, 'G': 290, 'H': 287, 'I': 284,
                'K': 224, 'L': 293, 'M': 281, 'N': 235,
                'P': 221, 'Q': 185, 'R': 244, 'S': 228,
                'T': 256, 'V': 315, 'W': 289, 'Y': 343
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def pka(
            self,
            kind='Taylor',
            standardize=True,
    ):
        """
        pKa: Negative of the logarithm of the acid dissociation
            constants for the COOH and NH2 groups (and, in
            some cases, other groups) in the molecule (at 25°C)

        References
        ----------
        Lundblad, R.L., & Macdonald, F. (Eds.). (2018). Handbook of Biochemistry
         and Molecular Biology (5th ed.). CRC Press. https://doi.org/10.1201/b21846

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Taylor':
            weight_dict = {
                'A': 2.33, 'C': 1.91, 'D': 1.95, 'E': 2.16,
                'F': 2.18, 'G': 2.34, 'H': 1.70, 'I': 2.26,
                'K': 2.15, 'L': 2.32, 'M': 2.16, 'N': 2.16,
                'P': 1.95, 'Q': 2.18, 'R': 2.03, 'S': 2.13,
                'T': 2.20, 'V': 2.27, 'W': 2.38, 'Y': 2.24
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def pkb(
            self,
            kind='Taylor',
            standardize=True,
    ):
        """
        pKa: Negative of the logarithm of the acid dissociation
            constants for the COOH and NH2 groups (and, in
            some cases, other groups) in the molecule (at 25°C)

        References
        ----------
        Lundblad, R.L., & Macdonald, F. (Eds.). (2018). Handbook of Biochemistry
         and Molecular Biology (5th ed.). CRC Press. https://doi.org/10.1201/b21846

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Taylor':
            weight_dict = {
                'A': 9.71, 'C': 10.28, 'D': 9.66, 'E': 9.58,
                'F': 9.09, 'G': 9.58, 'H': 9.09, 'I': 9.60,
                'K': 9.16, 'L': 9.58, 'M': 9.08, 'N': 8.73,
                'P': 10.47, 'Q': 9.00, 'R': 9.00, 'S': 9.05,
                'T': 8.96, 'V': 9.52, 'W': 9.34, 'Y': 9.04
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def hydrophilicity(
            self,
            kind='Hopp',
            standardize=True,
    ):
        """
        Hydrophilicity

        References
        ----------
        Hopp TP, Woods KR. Prediction of protein antigenic determinants from
        amino acid sequences. Proc Natl Acad Sci U S A. 1981 Jun;78(6):3824-8.
        doi: 10.1073/pnas.78.6.3824.

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Hopp':
            weight_dict = {
                'A': -0.5, 'C': -1.0, 'D': 3.0, 'E': 3.0,
                'F': -2.5, 'G': 0.0, 'H': -0.5, 'I': -1.8,
                'K': 3.0, 'L': -1.8, 'M': -1.3, 'N': 0.2,
                'P': 0.0, 'Q': 0.2, 'R': 3.0, 'S': 0.3,
                'T': -0.4, 'V': -1.5, 'W': -3.4, 'Y': -2.3
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def hydrophobicity(
            self,
            kind='Argos',
            standardize=True,
    ):
        """
        Hydrophobicity

        References
        ----------
        ARGOS, P., RAO, J.K.M. and HARGRAVE, P.A. (1982), Structural Prediction of Membrane-Bound
        Proteins. European Journal of Biochemistry, 128: 565-575.
        https://doi.org/10.1111/j.1432-1033.1982.tb07002.x

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Argos':
            weight_dict = {
                'A': 0.61, 'C': 1.07, 'D': 0.46, 'E': 0.47,
                'F': 2.02, 'G': 0.07, 'H': 0.61, 'I': 2.22,
                'K': 1.15, 'L': 1.53, 'M': 1.18, 'N': 0.06,
                'P': 1.95, 'Q': 0.00, 'R': 0.60, 'S': 0.05,
                'T': 0.05, 'V': 1.32, 'W': 2.65, 'Y': 1.88
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def fet(
            self,
            kind='Argos',
            standardize=True,
    ):
        """
        Fet

        References
        ----------
        ARGOS, P., RAO, J.K.M. and HARGRAVE, P.A. (1982), Structural Prediction of Membrane-Bound
        Proteins. European Journal of Biochemistry, 128: 565-575.
        https://doi.org/10.1111/j.1432-1033.1982.tb07002.x

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Argos':
            weight_dict = {
                'A': 1.24, 'C': 1.60, 'D': 0.71,
                'E': 0.65, 'F': 1.36, 'G': 1.24,
                'H': 1.01, 'I': 1.48, 'K': 0.00,
                'L': 1.36, 'M': 1.30, 'N': 0.77,
                'P': 0.89, 'Q': 0.65, 'R': 0.24,
                'S': 1.01, 'T': 0.95, 'V': 1.42,
                'W': 0.24, 'Y': 0.83
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict
        
    def hydration(
            self,
            kind='Argos',
            standardize=True,
    ):
        """
        Hydration

        References
        ----------
        ARGOS, P., RAO, J.K.M. and HARGRAVE, P.A. (1982), Structural Prediction of Membrane-Bound
        Proteins. European Journal of Biochemistry, 128: 565-575.
        https://doi.org/10.1111/j.1432-1033.1982.tb07002.x

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Argos':
            weight_dict = {
                'A': 2.00, 'C': 1.51, 'D': 0.01,
                'E': 0.13, 'F': 1.58, 'G': 2.07,
                'H': 0.12, 'I': 2.03, 'K': 0.23,
                'L': 2.06, 'M': 1.47, 'N': 0.20,
                'P': 0.93, 'Q': 0.25, 'R': 0.00,
                'S': 0.92, 'T': 0.94, 'V': 2.00,
                'W': 0.79, 'Y': 0.75
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def signal(
            self,
            kind='Argos',
            standardize=True,
    ):
        """
        Signal

        References
        ----------
        ARGOS, P., RAO, J.K.M. and HARGRAVE, P.A. (1982), Structural Prediction of Membrane-Bound
        Proteins. European Journal of Biochemistry, 128: 565-575.
        https://doi.org/10.1111/j.1432-1033.1982.tb07002.x

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Argos':
            weight_dict = {
                'A': 1.81, 'C': 1.89, 'D': 0.05,
                'E': 0.11, 'F': 1.96, 'G': 0.49,
                'H': 0.31, 'I': 1.45, 'K': 0.06,
                'L': 3.23, 'M': 2.67, 'N': 0.23,
                'P': 0.76, 'Q': 0.72, 'R': 0.20,
                'S': 0.97, 'T': 0.84, 'V': 1.08,
                'W': 0.77, 'Y': 0.39
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def volume(
            self,
            kind='Grantham',
            standardize=True,
    ):
        """
        Volume

        References
        ----------
        R. Grantham ,Amino Acid Difference Formula to Help Explain
        Protein Evolution.Science185,862-864(1974).DOI:10.1126/science.185.4154.862

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Grantham':
            weight_dict = {
                'A': 31, 'C': 55, 'D': 54, 'E': 83,
                'F': 132, 'G': 3, 'H': 96, 'I': 111,
                'K': 119, 'L': 111, 'M': 105, 'N': 56,
                'P': 32.5, 'Q': 85, 'R': 124, 'S': 32,
                'T': 61, 'V': 84, 'W': 170, 'Y': 136
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def polarity(
            self,
            kind='Grantham',
            standardize=True,
    ):
        """
        Polarity

        References
        ----------
        R. Grantham ,Amino Acid Difference Formula to Help Explain
        Protein Evolution.Science185,862-864(1974).DOI:10.1126/science.185.4154.862

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Grantham':
            weight_dict = {
                'A': 8.1, 'C': 5.5, 'D': 13.0, 'E': 12.3,
                'F': 5.2, 'G': 9.0, 'H': 10.4, 'I': 5.2,
                'K': 11.3, 'L': 4.9, 'M': 5.7, 'N': 11.6,
                'P': 8.0, 'Q': 10.5, 'R': 10.5, 'S': 9.2,
                'T': 8.6, 'V': 5.9, 'W': 5.4, 'Y': 6.2
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def composition(
            self,
            kind='Grantham',
            standardize=True,
    ):
        """
        Composition

        References
        ----------
        R. Grantham ,Amino Acid Difference Formula to Help Explain
        Protein Evolution.Science185,862-864(1974).DOI:10.1126/science.185.4154.862

        Parameters
        ----------
        kind
            method
        standardize
            if are values standalized

        Returns
        -------

        """
        if kind == 'Grantham':
            weight_dict = {
                'A': 0.00, 'C': 2.75, 'D': 1.38, 'E': 0.92,
                'F': 0.00, 'G': 0.74, 'H': 0.58, 'I': 0.00,
                'K': 0.33, 'L': 0.00, 'M': 0.00, 'N': 1.33,
                'P': 0.39, 'Q': 0.89, 'R': 0.65, 'S': 1.42,
                'T': 0.71, 'V': 0.00, 'W': 0.13, 'Y': 0.20
            }
            if standardize:
                return self.sstandardize.minmax2(weight_dict)
            else:
                return weight_dict

    def tobool(
            self,
            num_types,
            rank=None,
    ):
        """

        Parameters
        ----------
        num_types
            number of elements to be appended as bool
        rank
            rank of the 1st element to be appended

        Returns
        -------

        """
        bool_ = [0] * num_types
        if rank is not None:
            bool_[rank-1] = 1
            return bool_
        else:
            return bool_

    def get_(
            self,
            list_2d,
            window_aa_names,
    ):
        """
        ..  @summary:
            ---------
            Amino adcid descriptors are assembled.

        ..  @link:
            ------
            http://www.russelllab.org/aas/

        ..  @description:
            -------------
            > Charge, Hydrophobic, Aliphatic, Aromatic, Polar, Negative, Positive, Small
            are drawn by the substitution of amino acid property at @link.
            > Hydrophobicity, PI, Volume, Mass are values.

        ..  @since:
            -------
            It was introduced since v1.0

        :return: 2d dictionary
        :rtype: dict
        """
        start_time = time.time()
        list_2d_ = list_2d
        ### #/* block 1. */
        posi_ = self.positive()
        nega_ = self.negative()
        aliphatic_ = self.aliphatic()
        aromatic_ = self.aromatic()
        polar_ = self.polar()
        hydrophobic_ = self.hydrophobic()
        charged_ = self.charged()
        small_ = self.small()
        ### #/* block 2. */
        hydrophobicity = self.hydrophobicity()
        Hydrophilicity = self.hydrophilicity()
        volume = self.volume()
        aliphatic = self.weight()
        pI = self.pI()
        active = self.active()
        fet = self.fet()
        pka = self.pka()
        pkb = self.pkb()
        polarity = self.polarity()
        hydration = self.hydration()
        for i, aa_win_ids in enumerate(window_aa_names):
            for j in aa_win_ids:
                if j is None:
                    for e in self.tobool(num_types=3, rank=None):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=2, rank=None):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=2, rank=None):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=2, rank=None):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=3, rank=None):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=4, rank=None):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=4, rank=None):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=3, rank=None):
                        list_2d_[i].append(e)
                    list_2d_[i].append(0)
                    list_2d_[i].append(0)
                    list_2d_[i].append(0)
                    list_2d_[i].append(0)
                    list_2d_[i].append(0)
                    list_2d_[i].append(0)
                    list_2d_[i].append(0)
                    list_2d_[i].append(0)
                    list_2d_[i].append(0)
                    list_2d_[i].append(0)
                    list_2d_[i].append(0)
                else:
                    for e in self.tobool(num_types=3, rank=posi_[j]):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=2, rank=nega_[j]):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=2, rank=aliphatic_[j]):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=2, rank=aromatic_[j]):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=3, rank=polar_[j]):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=4, rank=hydrophobic_[j]):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=4, rank=charged_[j]):
                        list_2d_[i].append(e)
                    for e in self.tobool(num_types=3, rank=small_[j]):
                        list_2d_[i].append(e)
                    list_2d_[i].append(hydrophobicity[j])
                    list_2d_[i].append(Hydrophilicity[j])
                    list_2d_[i].append(volume[j])
                    list_2d_[i].append(aliphatic[j])
                    list_2d_[i].append(pI[j])
                    list_2d_[i].append(active[j])
                    list_2d_[i].append(fet[j])
                    list_2d_[i].append(pka[j])
                    list_2d_[i].append(pkb[j])
                    list_2d_[i].append(polarity[j])
                    list_2d_[i].append(hydration[j])
        # print(len(list_2d_[0]))
        # for i in list_2d_:
        #     if len(i) != 306:
        #         print(len(i))
        end_time = time.time()
        print('------> aa property {time}s.'.format(time=end_time - start_time))
        return list_2d_
    # print(p.get_(features, window_aa_names)[0])