This code supports the paper: 

Capra JA and Singh M. Predicting functionally important residues from
sequence conservation. Bioinformatics. 23(15):1875-82, 2007.  

If you use it, please cite the paper.

It is released under the GNU General Public License. Please see
COPYING for more details. It requires python and, if you wish to use
the von Neumann entropy method, numarray.  The other scoring methods
will work without numarray.  See the comments at the top of
score_conservation.py for more details.

Dependencies:
numarray (for vn_entropy only) -
http://sourceforge.net/project/showfiles.php?group_id=1369&release_id=223264

The program accepts alignments in the CLUSTAL and FASTA formats.  The
simplest usage is:

python score_conservation.py alignment.clustal

This will compute conservation scores for the alignment using the
Jensen-Shannon divergence with default settings and print out the
scores.

All methods described in the paper are implemented in
score_conservation.py. The relevant parameters can be specified with
command line options. For example, if you wanted to score an alignment
named alignment.clustal using Jensen-Shannon divergence, a window of
size 3 (on either side of the residue), and the swissprot background
distribution, and then save it in a file called alignment.scores, you
would type:

python score_conservation.py -s js_divergence -w 3 -d
      swissprot.distribution -o alignment.scores alignment.clustal

To get ouput in reference to a particular sequence in the alignment,
use the (-a) option.

Type: 

python score_conservation.py -h

for usage, the full list of options, and defaults. 

See the paper for a description of the methods. 

The .distribution files provide amino acid distributions. Feel free to
create your own. The matrix directory contains several BLOSUM matrices
and target frequencies. Again, any matrix file that follows the format
can be used.

