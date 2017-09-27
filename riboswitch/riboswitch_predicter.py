'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=too-few-public-methods
from __future__ import division

from collections import OrderedDict
from subprocess import Popen
import re
import sys
import tempfile

from synbiochem.utils.seq_utils import get_all_rev_trans, write_fasta
import pandas as pd


class RiboswitchPredicter(object):
    '''RiboswitchPredicter class.'''

    def __init__(self, pre_seqs, trunc_len, mutate_seq, post_seqs):
        self.__df = pd.DataFrame()
        self.__pre_seqs = pre_seqs if pre_seqs is not None else []
        self.__trunc_len = trunc_len
        self.__post_seqs = post_seqs if post_seqs is not None else []
        self.__df['variant'] = get_all_rev_trans(mutate_seq)
        self.__df['gc'] = _get_gc(self.__df['variant'])

    def get_data(self, temps=None):
        '''Gets data.'''
        for pre_idx in reversed(range(len(self.__pre_seqs))):
            pre_seq = ''.join(self.__pre_seqs[pre_idx:])

            for post_idx in range(len(self.__post_seqs) + 1):
                post_seq = ''.join(self.__post_seqs[:post_idx])
                seqs = self.__get_seqs(pre_seq, post_seq)
                trunc_seqs = self.__get_seqs(pre_seq[-self.__trunc_len:],
                                             post_seq)

                if temps is None:
                    temps = [37.0]

                for temp in temps:
                    orig_rnafold = _run_rnafold(seqs, temp=temp)
                    trunc_rnafold = _run_rnafold(trunc_seqs, temp=temp)

                    suf = '_'.join([str(val)
                                    for val in [len(pre_seq),
                                                len(self.__df['variant'][0]) +
                                                len(post_seq),
                                                temp]])

                    self.__df['dg_' + suf] = [val[1] for val in orig_rnafold]
                    self.__df['dg_trunc_' + suf] = \
                        [val[1] for val in trunc_rnafold]
                    self.__df['ddg_' + suf] = \
                        self.__df['dg_' + suf] - self.__df['dg_trunc_' + suf]

                    self.__df['struct_' + suf] = \
                        [val[0] for val in orig_rnafold]
                    self.__df['struct_trunc_' + suf] = \
                        [val[0] for val in trunc_rnafold]

        return self.__df

    def __get_seqs(self, pre_seq, post_seq):
        '''Get sequences.'''
        return [pre_seq + variant + post_seq
                for variant in self.__df['variant']]


def _run_rnafold(seqs, temp=37.0):
    '''Run RNAfold.'''
    fasta_filename = write_fasta(OrderedDict(zip(range(len(seqs)), seqs)))

    rnafold_file = tempfile.NamedTemporaryFile(prefix='rnafold_',
                                               suffix='.txt',
                                               delete=False)

    with open(fasta_filename) as fasta_file, \
            open(rnafold_file.name, 'w') as rnafold_file:

        # This calls RNAfold (kinda like from the command line):
        process = Popen(['RNAfold',
                         '--noPS',
                         '-T', str(temp)],
                        stdin=fasta_file,
                        stdout=rnafold_file)
        process.wait()
        rnafold_file.flush()

    # Read raw results file:
    return _read_rnafold_file(rnafold_file.name)


def _get_gc(seqs):
    '''Get GC content.'''
    return [(seq.count('G') + seq.count('C')) / len(seq) for seq in seqs]


def _read_rnafold_file(rnafold_filename):
    '''Read RNAfold file file.'''
    results = []

    with open(rnafold_filename) as rnafold_file:
        for line in rnafold_file:
            # Look to see if the line contains a number:
            mfe = re.findall(r'(?<=[^>])([+-]?\d+.\d+)', line)

            if mfe:
                results.append([line.split()[0], float(mfe[0])])

    return results


def main(args):
    '''main method.'''
    trunc_idx = -1

    for trunc_idx, arg in enumerate(args):
        try:
            int(arg)
            break
        except ValueError:
            continue

    pre_seqs = args[1:trunc_idx]
    trunc_len = int(args[trunc_idx])
    mutate_seq = args[trunc_idx + 1]
    post_seqs = args[trunc_idx + 2:]

    rib_pred = RiboswitchPredicter(pre_seqs, trunc_len, mutate_seq, post_seqs)
    df = rib_pred.get_data()
    print df
    df.to_csv(args[0] + '.csv', index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
