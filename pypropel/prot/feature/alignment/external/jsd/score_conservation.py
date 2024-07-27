#!/usr/bin/env python

################################################################################
# score_conservation.py - Copyright Tony Capra 2007 - Last Update: 03/09/11
#
# 03/09/11 - default window size = 3, as stated in help
# 06/21/09 - seq specific output now compatible with ConCavity
# 06/21/09 - numarray only included when vn_entropy is used
# 08/15/08 - added property_relative_entropy scoring method
# 08/15/08 - added equal sequence length check
# 01/07/08 - added z-score normalization option (-n)
# 01/07/08 - added seq. specific output option (-a)
# 11/30/07 - read_scoring_matrix now returns list rather than array.
# 11/30/07 - added window lambda command line option (-b)
# 07/05/07 - fixed gap penalty cutoff (<=) and error message
# 06/26/07 - fixed read_clustal_align bug
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# -----------------------------------------------------------------------------
#
# Dependencies: 
# numarray (for von Neumann entropy method) - 
# http://sourceforge.net/project/showfiles.php?group_id=1369&release_id=223264
# 
# 
# See usage() for usage.
#
# This program supports the paper: Capra JA and Singh M. 
# Predicting functionally important residues from sequence
# conservation. Bioinformatics. 23(15): 1875-1882, 2007.
# Please cite this paper if you use this code.
#
# It contains code for each method of scoring conservation that is
# evaluated: sum of pairs, weighted sum of pairs, Shannon entropy,
# Shannon entropy with property groupings (Mirny and Shakhnovich 95,
# Valdar and Thornton 01), relative entropy with property groupings
# (Williamson 95), von Neumann entropy (Caffrey et al 04), relative
# entropy (Samudrala and Wang 06), and Jensen-Shannon divergence
# (Capra and Singh 07).
#
# The code distributed with Mayrose et al 04 is used for Rate4Site. As
# of today it can be obtained from:
# http://www.tau.ac.il/~itaymay/cp/rate4site.html
# 
# 
# All scoring functions follow the same prototype:
#
# def score(col, sim_matrix, bg_disr, seq_weights, gap_penalty=1):
#
# - col: the column to be scored.
# 
# - sim_matrix: the similarity (scoring) matrix to be used. Not all 
# methods will use this parameter. 
# 
# - bg_distr: a list containing an amino acid probability distribution. Not
# all methods use this parameter. The default is the blosum62 background, but
# other distributions can be given. 
#
# - seq_weights: an array of floats that is used to weight the contribution
# of each seqeuence. If the len(seq_weights) != len(col), then every sequence 
# gets a weight of one.
#
# - gap_penalty: a binary variable: 0 for no gap penalty and 1
# for gap penalty. The default is to use a penalty. The gap penalty used is
# the score times the fraction of non-gap positions in the column.
#
#
# For a window score of any of above methods use the window_score method to 
# transform the individual column scores. 
#
################################################################################

import math, sys, getopt

# numarray imported below

PSEUDOCOUNT = .0000001

amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-']
iupac_alphabet = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "U", "V",
                  "W", "Y", "Z", "X", "*", "-"]

# dictionary to map from amino acid to its row/column in a similarity matrix
aa_to_index = {}
for i, aa in enumerate(amino_acids):
    aa_to_index[aa] = i


def usage():
    print("""\nUSAGE:\npython score_conservation.py [options] alignfile\n\t -alignfile must be in fasta or clustal format.\n\nOPTIONS:\n\t
    -a\treference sequence. Print scores in reference to a specific sequence (ignoring gaps). Default prints the entire column. [sequence name]\n\t
    -b\tlambda for window heuristic linear combination. Default=.5 [real in [0,1]]\n
    -d\tbackground distribution file, e.g., swissprot.distribution. Default=BLOSUM62 background [filename]\n\t
    -g\tgap cutoff. Do not score columns that contain more than gap cutoff fraction gaps. Default=.3 [real in [0, 1)]\n\t
    -h\thelp. Print this message.\n
    -l\tuse sequence weighting. Default=True [True|False]\n\t
    -m\tsimilarity matrix file, e.g., matrix/blosum62.bla or .qij. Default=identity matrix [filename]\n\t
    -n\tnormalize scores. Print the z-score (over the alignment) of each column raw score. Default=False\n\t
    -o\tname of output file. Default=output to screen [filename]\n\t
    -p\tuse gap penalty. Lower the score of columns that contain gaps. Default=True [True|False]\n\t
    -s\tconservation estimation method. \n\t\tOptions: shannon_entropy, property_entropy, property_relative_entropy, vn_entropy, relative_entropy, js_divergence, sum_of_pairs. Default=js_divergence\n\t
    -w\twindow size. Number of residues on either side included in the window. Default=3 [int]\n\t
    """)


################################################################################
# Frequency Count and Gap Penalty
################################################################################

def weighted_freq_count_pseudocount(col, seq_weights, pc_amount):
    """ Return the weighted frequency count for a column--with pseudocount."""

    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1.] * len(col)

    aa_num = 0
    freq_counts = len(amino_acids) * [pc_amount]  # in order defined by amino_acids

    for aa in amino_acids:
        for j in range(len(col)):
            if col[j] == aa:
                freq_counts[aa_num] += 1 * seq_weights[j]

        aa_num += 1

    for j in range(len(freq_counts)):
        freq_counts[j] = freq_counts[j] / (sum(seq_weights) + len(amino_acids) * pc_amount)

    return freq_counts


def weighted_gap_penalty(col, seq_weights):
    """ Calculate the simple gap penalty multiplier for the column. If the 
    sequences are weighted, the gaps, when penalized, are weighted 
    accordingly. """

    # if the weights do not match, use equal weight
    if len(seq_weights) != len(col):
        seq_weights = [1.] * len(col)

    gap_sum = 0.
    for i in range(len(col)):
        if col[i] == '-':
            gap_sum += seq_weights[i]

    return 1 - (gap_sum / sum(seq_weights))


def gap_percentage(col):
    """Return the percentage of gaps in col."""
    num_gaps = 0.

    for aa in col:
        if aa == '-': num_gaps += 1

    return num_gaps / len(col)


################################################################################
# Shannon Entropy
################################################################################

def shannon_entropy(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
    """Calculates the Shannon entropy of the column col. sim_matrix  and 
    bg_distr are ignored. If gap_penalty == 1, then gaps are penalized. The 
    entropy will be between zero and one because of its base. See p.13 of 
    Valdar 02 for details. The information score 1 - h is returned for the sake 
    of consistency with other scores."""

    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

    h = 0.
    for i in range(len(fc)):
        if fc[i] != 0:
            h += fc[i] * math.log(fc[i])

    #    h /= math.log(len(fc))
    h /= math.log(min(len(fc), len(col)))

    inf_score = 1 - (-1 * h)

    if gap_penalty == 1:
        return inf_score * weighted_gap_penalty(col, seq_weights)
    else:
        return inf_score


################################################################################
# Property Entropy
################################################################################

def property_entropy(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
    """Calculate the entropy of a column col relative to a partition of the 
    amino acids. Similar to Mirny '99. sim_matrix and bg_distr are ignored, but 
    could be used to define the sets. """

    # Mirny and Shakn. '99
    property_partition = [['A', 'V', 'L', 'I', 'M', 'C'], ['F', 'W', 'Y', 'H'], ['S', 'T', 'N', 'Q'], ['K', 'R'],
                          ['D', 'E'], ['G', 'P'], ['-']]

    # Williamson '95
    # property_partition = [['V','L', 'I','M'], ['F','W','Y'], ['S','T'], ['N','Q'], ['H','K','R'], ['D','E'], ['A','G'], ['P'], ['C'], ['-']]

    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

    # sum the aa frequencies to get the property frequencies
    prop_fc = [0.] * len(property_partition)
    for p in range(len(property_partition)):
        for aa in property_partition[p]:
            prop_fc[p] += fc[aa_to_index[aa]]

    h = 0.
    for i in range(len(prop_fc)):
        if prop_fc[i] != 0:
            h += prop_fc[i] * math.log(prop_fc[i])

    h /= math.log(min(len(property_partition), len(col)))

    inf_score = 1 - (-1 * h)

    if gap_penalty == 1:
        return inf_score * weighted_gap_penalty(col, seq_weights)
    else:
        return inf_score


################################################################################
# Property Relative Entropy
################################################################################

def property_relative_entropy(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
    """Calculate the relative entropy of a column col relative to a
    partition of the amino acids. Similar to Williamson '95. sim_matrix is
    ignored, but could be used to define the sets. See shannon_entropy()
    for more general info. """

    # Mirny and Shakn. '99
    # property_partition = [['A','V','L','I','M','C'], ['F','W','Y','H'], ['S','T','N','Q'], ['K','R'], ['D', 'E'], ['G', 'P'], ['-']]

    # Williamson '95
    property_partition = [['V', 'L', 'I', 'M'], ['F', 'W', 'Y'], ['S', 'T'], ['N', 'Q'], ['H', 'K', 'R'], ['D', 'E'],
                          ['A', 'G'], ['P'], ['C']]

    prop_bg_freq = []
    if len(bg_distr) == len(property_partition):
        prop_bg_freq = bg_distr
    else:
        prop_bg_freq = [0.248, 0.092, 0.114, 0.075, 0.132, 0.111, 0.161, 0.043, 0.024, 0.000]  # from BL62

    # fc = weighted_freq_count_ignore_gaps(col, seq_weights)
    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

    # sum the aa frequencies to get the property frequencies
    prop_fc = [0.] * len(property_partition)
    for p in range(len(property_partition)):
        for aa in property_partition[p]:
            prop_fc[p] += fc[aa_to_index[aa]]

    d = 0.
    for i in range(len(prop_fc)):
        if prop_fc[i] != 0 and prop_bg_freq[i] != 0:
            d += prop_fc[i] * math.log(prop_fc[i] / prop_bg_freq[i], 2)

    if gap_penalty == 1:
        return d * weighted_gap_penalty(col, seq_weights)
    else:
        return d


################################################################################
# von Neumann Entropy
################################################################################

def vn_entropy(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
    """ Calculate the von Neuman Entropy as described in Caffrey et al. 04.
    This code was adapted from the implementation found in the PFAAT project 
    available on SourceForge. bg_distr is ignored."""

    aa_counts = [0.] * 20
    for aa in col:
        if aa != '-': aa_counts[aa_to_index[aa]] += 1

    dm_size = 0
    dm_aas = []
    for i in range(len(aa_counts)):
        if aa_counts[i] != 0:
            dm_aas.append(i)
            dm_size += 1

    if dm_size == 0: return 0.0

    row_i = 0
    col_i = 0
    dm = zeros((dm_size, dm_size), Float32)
    for i in range(dm_size):
        row_i = dm_aas[i]
        for j in range(dm_size):
            col_i = dm_aas[j]
            dm[i][j] = aa_counts[row_i] * sim_matrix[row_i][col_i]

    ev = la.eigenvalues(dm).real

    temp = 0.
    for e in ev:
        temp += e

    if temp != 0:
        for i in range(len(ev)):
            ev[i] = ev[i] / temp

    vne = 0.0
    for e in ev:
        if e > (10 ** -10):
            vne -= e * math.log(e) / math.log(20)

    if gap_penalty == 1:
        # return (1-vne) * weighted_gap_penalty(col, seq_weights)
        return (1 - vne) * weighted_gap_penalty(col, [1.] * len(col))
    else:
        return 1 - vne


################################################################################
# Relative Entropy 
################################################################################

def relative_entropy(col, sim_matix, bg_distr, seq_weights, gap_penalty=1):
    """Calculate the relative entropy of the column distribution with a 
    background distribution specified in bg_distr. This is similar to the 
    approach proposed in Wang and Samudrala 06. sim_matrix is ignored."""

    distr = bg_distr[:]

    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

    # remove gap count
    if len(distr) == 20:
        new_fc = fc[:-1]
        s = sum(new_fc)
        for i in range(len(new_fc)):
            new_fc[i] = new_fc[i] / s
        fc = new_fc

    if len(fc) != len(distr): return -1

    d = 0.
    for i in range(len(fc)):
        if distr[i] != 0.0:
            d += fc[i] * math.log(fc[i] / distr[i])

    d /= math.log(len(fc))

    if gap_penalty == 1:
        return d * weighted_gap_penalty(col, seq_weights)
    else:
        return d


################################################################################
# Jensen-Shannon Divergence
################################################################################

def js_divergence(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
    """ Return the Jensen-Shannon Divergence for the column with the background
    distribution bg_distr. sim_matrix is ignored. JSD is the default method."""

    distr = bg_distr[:]

    fc = weighted_freq_count_pseudocount(col, seq_weights, PSEUDOCOUNT)

    # if background distrubtion lacks a gap count, remove fc gap count
    if len(distr) == 20:
        new_fc = fc[:-1]
        s = sum(new_fc)
        for i in range(len(new_fc)):
            new_fc[i] = new_fc[i] / s
        fc = new_fc

    if len(fc) != len(distr): return -1

    # make r distriubtion
    r = []
    for i in range(len(fc)):
        r.append(.5 * fc[i] + .5 * distr[i])

    d = 0.
    for i in range(len(fc)):
        if r[i] != 0.0:
            if fc[i] == 0.0:
                d += distr[i] * math.log(distr[i] / r[i], 2)
            elif distr[i] == 0.0:
                d += fc[i] * math.log(fc[i] / r[i], 2)
            else:
                d += fc[i] * math.log(fc[i] / r[i], 2) + distr[i] * math.log(distr[i] / r[i], 2)

    # d /= 2 * math.log(len(fc))
    d /= 2

    if gap_penalty == 1:
        return d * weighted_gap_penalty(col, seq_weights)
    else:
        return d


################################################################################
# Mutation Weighted Pairwise Match
################################################################################

def sum_of_pairs(col, sim_matrix, bg_distr, seq_weights, gap_penalty=1):
    """ Sum the similarity matrix values for all pairs in the column. 
    This method is similar to those proposed in Valdar 02. bg_distr is ignored."""

    sum = 0.
    max_sum = 0.

    for i in range(len(col)):
        for j in range(i):
            if col[i] != '-' and col[j] != '-':
                max_sum += seq_weights[i] * seq_weights[j]
                sum += seq_weights[i] * seq_weights[j] * sim_matrix[aa_to_index[col[i]]][aa_to_index[col[j]]]

    if max_sum != 0:
        sum /= max_sum
    else:
        sum = 0.

    if gap_penalty == 1:
        return sum * weighted_gap_penalty(col, seq_weights)
    else:
        return sum


################################################################################
# Window Score
################################################################################

def window_score(scores, window_len, lam=.5):
    """ This function takes a list of scores and a length and transforms them 
    so that each position is a weighted average of the surrounding positions. 
    Positions with scores less than zero are not changed and are ignored in the 
    calculation. Here window_len is interpreted to mean window_len residues on 
    either side of the current residue. """

    w_scores = scores[:]

    for i in range(window_len, len(scores) - window_len):
        if scores[i] < 0:
            continue

        sum = 0.
        num_terms = 0.
        for j in range(i - window_len, i + window_len + 1):
            if i != j and scores[j] >= 0:
                num_terms += 1
                sum += scores[j]

        if num_terms > 0:
            w_scores[i] = (1 - lam) * (sum / num_terms) + lam * scores[i]

    return w_scores


def calc_z_scores(scores, score_cutoff):
    """Calculates the z-scores for a set of scores. Scores below
    score_cutoff are not included."""

    average = 0.
    std_dev = 0.
    z_scores = []
    num_scores = 0

    for s in scores:
        if s > score_cutoff:
            average += s
            num_scores += 1
    if num_scores != 0:
        average /= num_scores

    for s in scores:
        if s > score_cutoff:
            std_dev += ((s - average) ** 2) / num_scores
    std_dev = math.sqrt(std_dev)

    for s in scores:
        if s > score_cutoff and std_dev != 0:
            z_scores.append((s - average) / std_dev)
        else:
            z_scores.append(-1000.0)

    return z_scores


################################################################################
################################################################################
################################################################################
#  END CONSERVATION SCORES
################################################################################
################################################################################
################################################################################


def read_scoring_matrix(sm_file):
    """ Read in a scoring matrix from a file, e.g., blosum80.bla, and return it
    as an array. """
    aa_index = 0
    first_line = 1
    row = []
    list_sm = []  # hold the matrix in list form

    try:
        matrix_file = open(sm_file, 'r')

        for line in matrix_file:

            if line[0] != '#' and first_line:
                first_line = 0
                if len(amino_acids) == 0:
                    for c in line.split():
                        aa_to_index[string.lower(c)] = aa_index
                        amino_acids.append(string.lower(c))
                        aa_index += 1

            elif line[0] != '#' and first_line == 0:
                if len(line) > 1:
                    row = line.split()
                    list_sm.append(row)

    except IOError:
        print("Could not load similarity matrix: %s. Using identity matrix..." % sm_file)
        return identity(20)

    # if matrix is stored in lower tri form, copy to upper
    if len(list_sm[0]) < 20:
        for i in range(0, 19):
            for j in range(i + 1, 20):
                list_sm[i].append(list_sm[j][i])

    for i in range(len(list_sm)):
        for j in range(len(list_sm[i])):
            list_sm[i][j] = float(list_sm[i][j])

    return list_sm
    # sim_matrix = array(list_sm, type=Float32)
    # return sim_matrix


def calculate_sequence_weights(msa):
    """ Calculate the sequence weights using the Henikoff '94 method
    for the given msa. """

    seq_weights = [0.] * len(msa)
    for i in range(len(msa[0])):
        freq_counts = [0] * len(amino_acids)

        col = []
        for j in range(len(msa)):
            if msa[j][i] != '-':  # ignore gaps
                freq_counts[aa_to_index[msa[j][i]]] += 1

        num_observed_types = 0
        for j in range(len(freq_counts)):
            if freq_counts[j] > 0: num_observed_types += 1

        for j in range(len(msa)):
            d = freq_counts[aa_to_index[msa[j][i]]] * num_observed_types
            if d > 0:
                seq_weights[j] += 1. / d

    for w in range(len(seq_weights)):
        seq_weights[w] /= len(msa[0])

    return seq_weights


def load_sequence_weights(fname):
    """Read in a sequence weight file f and create sequence weight list. 
    The weights are in the same order as the sequences each on a new line. """
    seq_weights = []

    try:
        f = open(fname)

        for line in f:
            l = line.split()
            if line[0] != '#' and len(l) == 2:
                seq_weights.append(float(l[1]))

    except IOError:
        pass
    # print "No sequence weights. Can't find: ", fname

    return seq_weights


def get_column(col_num, alignment):
    """Return the col_num column of alignment as a list."""
    col = []
    for seq in alignment:
        if col_num < len(seq): col.append(seq[col_num])

    return col


def get_distribution_from_file(fname):
    """ Read an amino acid distribution from a file. The probabilities should
    be on a single line separated by whitespace in alphabetical order as in 
    amino_acids above. # is the comment character."""

    distribution = []
    try:
        f = open(fname)
        for line in f:
            if line[0] == '#': continue
            line = line[:-1]
            distribution = line.split()
            distribution = map(float, distribution)


    except IOError:
        print ("Using default (BLOSUM62) background.")
        return []

    # use a range to be flexible about round off
    if .997 > sum(distribution) or sum(distribution) > 1.003:
        print("Distribution does not sum to 1. Using default (BLOSUM62) background.")
        print(sum(distribution))
        return []

    return distribution


def read_fasta_alignment(filename):
    """ Read in the alignment stored in the FASTA file, filename. Return two
    lists: the identifiers and sequences. """

    f = open(filename)

    names = []
    alignment = []
    cur_seq = ''

    for line in f:
        line = line[:-1]
        if len(line) == 0: continue

        if line[0] == ';': continue
        if line[0] == '>':
            names.append(line[1:].replace('\r', ''))

            if cur_seq != '':
                cur_seq = cur_seq.upper()
                for i, aa in enumerate(cur_seq):
                    if aa not in iupac_alphabet:
                        cur_seq = cur_seq.replace(aa, '-')
                alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))
                cur_seq = ''
        elif line[0] in iupac_alphabet:
            cur_seq += line.replace('\r', '')

    # add the last sequence
    cur_seq = cur_seq.upper()
    for i, aa in enumerate(cur_seq):
        if aa not in iupac_alphabet:
            cur_seq = cur_seq.replace(aa, '-')
    alignment.append(cur_seq.replace('B', 'D').replace('Z', 'Q').replace('X', '-'))

    return names, alignment


def read_clustal_alignment(filename):
    """ Read in the alignment stored in the CLUSTAL file, filename. Return
    two lists: the names and sequences. """

    names = []
    alignment = []

    f = open(filename)

    for line in f:
        line = line[:-1]
        if len(line) == 0: continue
        if '*' in line: continue

        if 'CLUSTAL' in line: continue

        t = line.split()

        if len(t) == 2 and t[1][0] in iupac_alphabet:
            if t[0] not in names:
                names.append(t[0])
                alignment.append(t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X', '-').replace('\r', ''))
            else:
                alignment[names.index(t[0])] += t[1].upper().replace('B', 'D').replace('Z', 'Q').replace('X',
                                                                                                         '-').replace(
                    '\r', '')

    return names, alignment


################################################################################
# Begin execution
################################################################################

# BLOSUM62 background distribution
blosum_background_distr = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024,
                           0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]

# set defaults
window_size = 3  # 0 = no window
win_lam = .5  # for window method linear combination
outfile_name = ""
import os
s_matrix_file = os.path.dirname(os.getcwd()) + '/alignment/external/jsd/' + "matrix/blosum62.bla"
bg_distribution = blosum_background_distr[:]
scoring_function = js_divergence
use_seq_weights = True
background_name = 'blosum62'
gap_cutoff = .3
use_gap_penalty = 1
seq_specific_output = 0  # name of sequence if True
normalize_scores = False

# parse options and args -- see usage()
if len(sys.argv) < 2:
    usage()
    sys.exit(2)

try:
    opts, args = getopt.getopt(sys.argv[1:], "hl:d:m:o:s:w:g:p:b:a:n")
except getopt.GetoptError:
    usage()
    sys.exit(1)

if len(args) < 1:
    usage()
    sys.exit(1)

for opt, arg in opts:
    if opt == "-h":
        usage()
        sys.exit()
    if opt == "-o":
        outfile_name = arg
    elif opt == "-l":
        if 'false' in arg.lower():
            use_seq_weights = False
    elif opt == "-p":
        if 'false' in arg.lower():
            use_gap_penalty = 0
    elif opt == "-m":
        s_matrix_file = arg
    elif opt == "-d":
        d = get_distribution_from_file(arg)
        if d != []:
            bg_distribution = d
            background_name = arg
    elif opt == "-w":
        try:
            window_size = int(arg)
        except ValueError:
            print("ERROR: Window size must be an integer. Using window_size 3...")
            window_size = 3
    elif opt == "-b":
        try:
            win_lam = float(arg)
            if not (0. <= win_lam <= 1.): raise ValueError
        except ValueError:
            print("ERROR: Window lambda must be a real in [0,1]. Using lambda = .5...")
            win_lam = .5
    elif opt == "-g":
        try:
            gap_cutoff = float(arg)
            if not (0. <= gap_cutoff < 1.): raise ValueError
        except ValueError:
            print("ERROR: Gap cutoff must be a real in [0,1). Using a gap cutoff of .3...")
            gap_cutoff = .3
    elif opt == '-a':
        seq_specific_output = arg
    elif opt == '-n':
        normalize_scores = True
    elif opt == '-s':
        if arg == 'shannon_entropy':
            scoring_function = shannon_entropy
        elif arg == 'property_entropy':
            scoring_function = property_entropy
        elif arg == 'property_relative_entropy':
            scoring_function = property_relative_entropy
        elif arg == 'vn_entropy':
            scoring_function = vn_entropy; from numarray import *; import numarray.linear_algebra as la

        elif arg == 'relative_entropy':
            scoring_function = relative_entropy
        elif arg == 'js_divergence':
            scoring_function = js_divergence
        elif arg == 'sum_of_pairs':
            scoring_function = sum_of_pairs
        else:
            print("%s is not a valid scoring method. Using %s.\n" % (arg, scoring_function.__name__))

align_file = args[0]
align_suffix = align_file.split('.')[-1]

s_matrix = read_scoring_matrix(s_matrix_file)

names = []
alignment = []
seq_weights = []

try:
    names, alignment = read_clustal_alignment(align_file)
    if names == []:
        names, alignment = read_fasta_alignment(align_file)
except IOError:
    print ("Could not find %s. Exiting..." % align_file)
    sys.exit(1)

if len(alignment) != len(names) or alignment == []:
    print("Unable to parse alignment.\n")
    sys.exit(1)

seq_len = len(alignment[0])
for i, seq in enumerate(alignment):
    if len(seq) != seq_len:
        print("ERROR: Sequences of different lengths: %s (%d) != %s (%d).\n" % (names[0], seq_len, names[i], len(seq)))
        sys.exit(1)

if use_seq_weights:
    seq_weights = load_sequence_weights(align_file.replace('.%s' % align_suffix, '.weights'))
    if seq_weights == []:
        seq_weights = calculate_sequence_weights(alignment)

if len(seq_weights) != len(alignment): seq_weights = [1.] * len(alignment)

# handle print of output relative to specific sequence
ref_seq_num = None
if seq_specific_output and seq_specific_output not in names:
    print("Sequence %s not found in alignment. Using default output format...\n" % seq_specific_output)
    seq_specific_output = 0
elif seq_specific_output in names:
    ref_seq_num = names.index(seq_specific_output)

# calculate scores
scores = []
for i in range(len(alignment[0])):
    col = get_column(i, alignment)

    if len(col) == len(alignment):
        if gap_percentage(col) <= gap_cutoff:
            scores.append(scoring_function(col, s_matrix, bg_distribution, seq_weights, use_gap_penalty))
        else:
            scores.append(-1000.)

if window_size > 0:
    scores = window_score(scores, window_size, win_lam)

if normalize_scores:
    scores = calc_z_scores(scores, -999)

# print to file/stdout
try:
    if outfile_name != "":
        outfile = open(outfile_name, 'w')
        outfile.write(
            "# %s -- %s - window_size: %d - window lambda: %.2f - background: %s - seq. weighting: %s - gap penalty: %d - normalized: %s\n" % (
            align_file, scoring_function.__name__, window_size, win_lam, background_name, use_seq_weights,
            use_gap_penalty, normalize_scores))
        if seq_specific_output:
            outfile.write("# reference sequence: %s\n" % seq_specific_output)
            outfile.write("# align_column_number\tamino acid\tscore\n")
        else:
            outfile.write("# align_column_number\tscore\tcolumn\n")
    else:
        print(
            "# %s -- %s - window_size: %d - background: %s - seq. weighting: %s - gap penalty: %d - normalized: %s" % (
            align_file, scoring_function.__name__, window_size, background_name, use_seq_weights, use_gap_penalty,
            normalize_scores))
        if seq_specific_output:
            print("# reference sequence: %s" % seq_specific_output)
            print("# align_column_number\tamino acid\tscore\n")
        else:
            print("# align_column_number\tscore\tcolumn\n")

except IOError:
    print("Could not open %s for output. Printing results to standard out..." % outfile_name)
    outfile_name = ""

for i, score in enumerate(scores):
    if seq_specific_output:
        cur_aa = get_column(i, alignment)[ref_seq_num]
        if cur_aa == '-': continue
        if outfile_name == "":
            print("%d\t%s\t%.5f" % (i, cur_aa, score))
        else:
            outfile.write("%d\t%s\t%5f\n" % (i, cur_aa, score))
    else:
        if outfile_name == "":
            print("%d\t%.5f\t%s" % (i, score, "".join(get_column(i, alignment))))
        else:
            outfile.write("%d\t%5f\t%s\n" % (i, score, "".join(get_column(i, alignment))))
