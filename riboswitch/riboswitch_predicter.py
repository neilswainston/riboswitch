'''
synbiochemdev (c) University of Manchester 2015

synbiochemdev is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
import sys

from riboswitch.vienna import get_minimum_free_energy
from synbiochem.utils import seq_utils


def analyse(pre_seq, trunc_len, tag, post_seq, outfile_name):
    '''Analyse.'''
    seqs = []
    pre_trunc_seq = pre_seq[-trunc_len:]

    for rev_trans in seq_utils.get_all_rev_trans(tag):
        seqs.extend([pre_seq + rev_trans + post_seq,
                     pre_trunc_seq + rev_trans + post_seq])

    mfes = get_minimum_free_energy(seqs)

    outfile = open(outfile_name, 'w')

    for i in xrange(0, len(seqs), 2):
        outfile.write('\t'.join([seqs[i], seqs[i + 1],
                                 str(mfes[i]), str(mfes[i + 1]),
                                 str(mfes[i] - mfes[i + 1])]) + '\n')

    outfile.close()


def main(args):
    '''main method'''
    pre_seq = args[0].strip()
    trunc_len = int(args[1].strip())
    tag = args[2].strip()
    post_seq = args[3].strip()
    outfile_name = args[4].strip()
    analyse(pre_seq, trunc_len, tag, post_seq, outfile_name)


if __name__ == '__main__':
    main(sys.argv[1:])
