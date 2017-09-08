'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from subprocess import Popen
import os
import re
import tempfile

from synbiochem.utils import seq_utils


def write_csv_to_fasta(csv_filename, fasta_filename):
    '''Writes a csv file of name, sequence to a fasta file.'''

    # Open csv file in read mode, fasta file in write mode:
    with open(csv_filename, 'r') as csv_file, \
            open(fasta_filename, 'w') as fasta_file:
        # Read csv file line by line...
        for line in csv_file:
            # Split line into "tokens", separated by commas:
            tokens = line.split(',')

            # Write header to fasta:
            fasta_file.write('>' + tokens[0] + '\n')

            # Write sequence to fasta:
            fasta_file.write(tokens[1])


def run_rnafold(fasta_filename):
    '''Run RNAfold from fasta file.'''

    # Generate raw results file:
    raw_filename = 'raw.txt'

    with open(fasta_filename) as fasta_file, \
            open(raw_filename, 'w') as raw_file:

        # This calls RNAfold (kinda like from the command line):
        process = Popen(['RNAfold', '--noPS'],
                        stdin=fasta_file,
                        stdout=raw_file)
        process.wait()
        raw_file.flush()

    # Read raw results file:
    results = {}

    with open(raw_filename) as raw_file:
        for line in raw_file:
            if line.startswith('>'):
                # If this is a fasta header line, store the name:
                name = line[1:].strip()
            else:
                # Look to see if the line contains a number:
                numbers = re.findall(r'[+-]?\d+.\d+', line)

                if numbers:
                    # Store name and number in results:
                    results[name] = float(numbers[0])

    return results


def get_minimum_free_energy(seqs):
    '''Returns minimum free energy of supplied DNA / RNA sequences.'''
    input_filename = seq_utils.write_fasta({str(idx): seq
                                            for idx, seq in enumerate(seqs)})

    with open(tempfile.NamedTemporaryFile().name, 'w') as output_file:
        proc = Popen('RNAfold',
                     stdin=open(input_filename),
                     stdout=output_file)

        proc.wait()

        _cleanup(os.getcwd(), 'Seq\\d+_ss.ps')

        mfes = []
        pattern = re.compile(r'[+-]?\d+\.\d+')

        with open(output_file.name) as out_file:
            for line in out_file.readlines():
                src = pattern.search(line)

                if src:
                    mfes.append(float(src.group()))

        return mfes


def write_results_to_csv(results, csv_filename):
    '''Writes the results to a csv file.'''
    with open(csv_filename, 'w') as csv_file:
        for name, value in results.iteritems():
            csv_file.write(name + ',' + str(value) + '\n')


def _cleanup(drctry, pattern):
    '''Deletes files in directory matching pattern.'''
    for filename in os.listdir(drctry):
        if re.search(pattern, filename):
            os.remove(os.path.join(drctry, filename))


def main():
    '''main method to start the program.'''
    fasta_filename = 'all_seq.fasta.txt'

    # Write fasta file:
    write_csv_to_fasta('all_seq.csv', fasta_filename)

    # Runs the RNAfold program with input from fasta file:
    results = run_rnafold(fasta_filename)

    # Writes the results to a csv file:
    write_results_to_csv(results, 'results.csv')


if __name__ == '__main__':
    main()
