__author__ = "Jianfeng Sun"
__version__ = "v1.0"
__copyright__ = "Copyright 2024"
__license__ = "GPL v3.0"
__email__ = "jianfeng.sunmt@gmail.com"
__maintainer__ = "Jianfeng Sun"

import numpy as np
import pandas as pd
from pypropel.prot.sequence.Symbol import Symbol


class Single:

    def __init__(
            self,
            msa,
    ):
        self.msa = msa
        self.aa_alphabet = Symbol().single(gap=True, universal=False)
        self.msa_row = len(self.msa)
        self.msa_col = len(self.msa[0])

    def columns(self, ):
        """
        Frequency of 20 amino acids and 1 gap in each column of msa.

        qi represents one of 20 amino acids and 1 gap in a column of msa.
           sum_t(qi) is total amount of qi in column t of msa. The total
           amount of rows in the msa is sum(rows). Frequencies matrix is
           calculated by sum(qi)/sum(rows).

        Examples
        --------
           1atzA: QPLDVILLLDGSSSFPASYFDEMKSFAKAFISKANIGPRLTQVSVLQYGSITTIDVPWNVVPEKAHLLSLVDVMQ
           return 75x21 matrix.

        Returns
        -------
        2d array : numpy.ndarray
            row: 75; col: 21

        """
        A = [0] * self.msa_col
        C = [0] * self.msa_col
        D = [0] * self.msa_col
        E = [0] * self.msa_col
        F = [0] * self.msa_col
        G = [0] * self.msa_col
        H = [0] * self.msa_col
        I = [0] * self.msa_col
        K = [0] * self.msa_col
        L = [0] * self.msa_col
        M = [0] * self.msa_col
        N = [0] * self.msa_col
        P = [0] * self.msa_col
        Q = [0] * self.msa_col
        R = [0] * self.msa_col
        S = [0] * self.msa_col
        T = [0] * self.msa_col
        V = [0] * self.msa_col
        W = [0] * self.msa_col
        Y = [0] * self.msa_col
        omit = [0] * self.msa_col
        # print(self.msa)
        for homolog in self.msa:
            # print(homolog)
            for alignment_pos, base in enumerate(homolog):
                if base == 'A':
                    A[alignment_pos] += 1
                elif base == 'C':
                    C[alignment_pos] += 1
                elif base == 'D':
                    D[alignment_pos] += 1
                elif base == 'E':
                    E[alignment_pos] += 1
                elif base == 'F':
                    F[alignment_pos] += 1
                elif base == 'G':
                    G[alignment_pos] += 1
                elif base == 'H':
                    H[alignment_pos] += 1
                elif base == 'I':
                    I[alignment_pos] += 1
                elif base == 'K':
                    K[alignment_pos] += 1
                elif base == 'L':
                    L[alignment_pos] += 1
                elif base == 'M':
                    M[alignment_pos] += 1
                elif base == 'N':
                    N[alignment_pos] += 1
                elif base == 'P':
                    P[alignment_pos] += 1
                elif base == 'Q':
                    Q[alignment_pos] += 1
                elif base == 'R':
                    R[alignment_pos] += 1
                elif base == 'S':
                    S[alignment_pos] += 1
                elif base == 'T':
                    T[alignment_pos] += 1
                elif base == 'V':
                    V[alignment_pos] += 1
                elif base == 'W':
                    W[alignment_pos] += 1
                elif base == 'Y':
                    Y[alignment_pos] += 1
                elif base == '-':
                    omit[alignment_pos] += 1
        count_array = np.array([A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, omit])
        # count_array = [A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, omit]
        # freq_array = count_array
        freq_array = count_array / self.msa_row
        freq = {}
        for i, aa in enumerate(self.aa_alphabet):
            freq[aa] = freq_array[i]
        return freq, np.array(omit) / self.msa_row

    def columnsByPandas(self, ):
        msa_sp = []
        for homolog in self.msa:
            msa_sp.append(list(homolog))
        msa_sp_df = pd.DataFrame(msa_sp)
        freq = {}
        for aa in self.aa_alphabet:
            freq[aa] = []
        for alignment_pos in msa_sp_df.columns:
            base_count = msa_sp_df[alignment_pos].value_counts().to_dict()
            for i, aa in enumerate(self.aa_alphabet):
                if aa in base_count.keys():
                    # freq[aa].append(base_count[aa] / self.msa_row)
                    freq[aa].append(base_count[aa])
                else:
                    freq[aa].append(0)
        return freq

    def alignment(self, ):
        """
        Frequency of 20 amino acids and 1 gap in a whole msa.

        qi represents one of 21 amino acids. sum(qi) is total amount of qi in a whole msa.
           The total amount of all amino acids in a whole msa is sum(qi). Frequencies are
           calculated by sum(qi)/sum(all).

        Examples
        --------
        1atzA: QPLDVILLLDGSSSFPASYFDEMKSFAKAFISKANIGPRLTQVSVLQYGSITTIDVPWNVVPEKAHLLSLVDVMQ
           return 1x21 matrix.

        Returns
        -------
            1d array: numpy.ndarray
             row: 1; col: 21
        """
        total_num_MSA = self.msa_row * self.msa_col
        A = 0
        C = 0
        D = 0
        E = 0
        F = 0
        G = 0
        H = 0
        I = 0
        K = 0
        L = 0
        M = 0
        N = 0
        P = 0
        Q = 0
        R = 0
        S = 0
        T = 0
        V = 0
        W = 0
        Y = 0
        omit = 0
        for i in range(self.msa_row):
            for j in range(self.msa_col):
                if self.msa[i][j] == 'A':
                    A += 1
                elif self.msa[i][j] == 'C':
                    C += 1
                elif self.msa[i][j] == 'D':
                    D += 1
                elif self.msa[i][j] == 'E':
                    E += 1
                elif self.msa[i][j] == 'F':
                    F += 1
                elif self.msa[i][j] == 'G':
                    G += 1
                elif self.msa[i][j] == 'H':
                    H += 1
                elif self.msa[i][j] == 'I':
                    I += 1
                elif self.msa[i][j] == 'K':
                    K += 1
                elif self.msa[i][j] == 'L':
                    L += 1
                elif self.msa[i][j] == 'M':
                    M += 1
                elif self.msa[i][j] == 'N':
                    N += 1
                elif self.msa[i][j] == 'P':
                    P += 1
                elif self.msa[i][j] == 'Q':
                    Q += 1
                elif self.msa[i][j] == 'R':
                    R += 1
                elif self.msa[i][j] == 'S':
                    S += 1
                elif self.msa[i][j] == 'T':
                    T += 1
                elif self.msa[i][j] == 'V':
                    V += 1
                elif self.msa[i][j] == 'W':
                    W += 1
                elif self.msa[i][j] == 'Y':
                    Y += 1
                elif self.msa[i][j] == '-':
                    omit += 1
        freq_array = np.array([
            A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y, omit
        ]) / total_num_MSA
        freq = {}
        for i, aa in enumerate(self.aa_alphabet):
            freq[aa] = freq_array[i]
        return freq