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
from sklearn.linear_model import LinearRegression
from sklearn.metrics.regression import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.preprocessing.imputation import Imputer

import matplotlib.pyplot as plt
import pandas as pd


def clean(df, strategy='median'):
    '''Cleans DataFrame.'''
    df.info()
    imputer = Imputer(strategy=strategy)
    object_df = df.select_dtypes(include=['object'])
    float_df = df.select_dtypes(include=['float64'])
    imputer.fit(float_df)
    float_df = pd.DataFrame(imputer.transform(float_df),
                            columns=float_df.columns)

    return pd.concat([object_df, float_df], axis=1)


def _plot_scatter(df):
    '''Make scatter plot.'''
    attributes = ['dg_105_18_37.0',
                  'dg_105_60_37.0',
                  'dg_105_844_37.0',
                  'dg_159_18_37.0',
                  'dg_159_60_37.0',
                  'dg_159_844_37.0']

    scatter_matrix(df[attributes], figsize=(8, 6))
    plt.show()


def _linear_regression(df):
    '''Perform linear regression.'''
    train, test = train_test_split(df, test_size=0.2)

    attributes = ['ddg_105_18_37.0',
                  'ddg_105_60_37.0',
                  'ddg_105_844_37.0',
                  'ddg_159_18_37.0',
                  'ddg_159_60_37.0',
                  'ddg_159_844_37.0']

    print train[attributes]

    lin_reg = LinearRegression()
    lin_reg.fit([[val] for val in train['ddg_105_18_37.0']],
                [[val] for val in train['max_min_ratio']])
    predictions = lin_reg.predict([[val] for val in test['ddg_105_18_37.0']])
    print predictions
    print test['max_min_ratio']
    print mean_squared_error([[val] for val in test['max_min_ratio']],
                             predictions)


def main(args):
    '''main method.'''
    df = pd.read_csv(args[0])
    results_df = clean(pd.read_csv(args[1]))
    results_df['max_min_ratio'] = \
        results_df['max_mean'] / results_df['min_mean']

    df = pd.merge(df, results_df, on='variant')
    df.to_csv('out.csv', index=False)

    # print df.describe()
    # df.hist(bins=10, figsize=(20, 15))
    # plt.show()

    # _plot_scatter(df)
    _linear_regression(df)


if __name__ == '__main__':
    main(sys.argv[1:])
