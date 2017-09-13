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

    def __init__(self, pre_seq, trunc_len, mutate_seq, post_seq):
        self.__df = pd.DataFrame()
        self.__pre_seq = pre_seq
        self.__trunc_len = trunc_len
        self.__post_seq = post_seq if post_seq is not None else ''
        self.__df['variant'] = get_all_rev_trans(mutate_seq)

    def get_data(self, temps=None):
        '''Gets data.'''
        seqs = self.__get_seqs(self.__pre_seq)
        trunc_seqs = self.__get_seqs(self.__pre_seq[-self.__trunc_len:])

        if temps is None:
            temps = [30.0, 37.0]

        for temp in temps:
            suf = str(temp)
            self.__df['dg_' + suf] = \
                _get_mfes(seqs, temp=temp)
            self.__df['dg_trunc_' + suf] = \
                _get_mfes(trunc_seqs, temp=temp)
            self.__df['ddg_' + suf] = \
                self.__df['dg_' + suf] - self.__df['dg_trunc_' + suf]

        self.__df['gc'] = _get_gc(seqs)
        self.__df['gc_trunc'] = _get_gc(trunc_seqs)

        return self.__df

    def __get_seqs(self, pre_seq):
        '''Get sequences.'''
        return [pre_seq + variant + self.__post_seq
                for variant in self.__df['variant']]


def _get_mfes(seqs, temp=37.0):
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
    return _read_rnafold_file(rnafold_file.name).values()


def _get_gc(seqs):
    '''Get GC content.'''
    return [(seq.count('G') + seq.count('C')) / len(seq) for seq in seqs]


def _read_rnafold_file(rnafold_filename):
    '''Read RNAfold file file.'''
    results = OrderedDict()

    with open(rnafold_filename) as rnafold_file:
        for line in rnafold_file:
            if line.startswith('>'):
                # If this is a fasta header line, store the name:
                name = line[1:].strip()
            else:
                # Look to see if the line contains a number:
                numbers = re.findall(r'[+-]?\d+.\d+', line)

                if numbers:
                    # Store name and number in results:
                    results[int(name)] = float(numbers[0])

    return results


def main(args):
    '''main method.'''
    pre_seq = args[0]
    trunc_len = int(args[1])
    mutate_seq = args[2]
    post_seq = args[3] if len(args) == 3 else None

    rib_pred = RiboswitchPredicter(pre_seq, trunc_len, mutate_seq, post_seq)
    df = rib_pred.get_data()
    print df
    df.to_csv(args[4] + '.csv', index=False)

    # Normalise:
    df_norm = df.ix[:, 1:]
    df_norm = (df_norm - df_norm.mean()) / df_norm.std()
    df_norm.insert(0, df.ix[:, 0].name, df.ix[:, 0].values)
    print df_norm
    df_norm.to_csv(args[4] + '_norm.csv', index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
