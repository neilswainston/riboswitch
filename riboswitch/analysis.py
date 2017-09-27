'''
synbiochem (c) University of Manchester 2017

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
# pylint: disable=ungrouped-imports
import sys

from pandas.plotting import scatter_matrix

import matplotlib.pyplot as plt
import pandas as pd


def main(args):
    '''main method.'''
    df = pd.read_csv(args[0])
    df.info()
    # print df.describe()
    # df.hist(bins=10, figsize=(20, 15))
    # plt.show()

    attributes = ['dg_105_18_37.0',
                  'dg_105_60_37.0',
                  'dg_105_844_37.0',
                  'dg_159_18_37.0',
                  'dg_159_60_37.0',
                  'dg_159_844_37.0']

    scatter_matrix(df[attributes], figsize=(8, 6))
    plt.show()


if __name__ == '__main__':
    main(sys.argv[1:])
