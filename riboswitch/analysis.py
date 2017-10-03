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
from sklearn.feature_selection import RFE
from sklearn.linear_model import LinearRegression, RandomizedLasso
from sklearn.metrics.regression import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing.imputation import Imputer

import matplotlib.pyplot as plt
import pandas as pd


def clean(df, strategy='median'):
    '''Cleans DataFrame.'''
    imputer = Imputer(strategy=strategy)
    object_df = df.select_dtypes(include=['object'])
    float_df = df.select_dtypes(include=['float64'])
    imputer.fit(float_df)
    float_df = pd.DataFrame(imputer.transform(float_df),
                            columns=float_df.columns)

    return pd.concat([object_df, float_df], axis=1)


def _plot_scatter(df, out='scatter.png'):
    '''Make scatter plot.'''
    attributes = ['gc',
                  'ddg_105_18_37.0',
                  'ddg_105_60_37.0',
                  'ddg_105_844_37.0',
                  'ddg_159_18_37.0',
                  'ddg_159_60_37.0',
                  'ddg_159_844_37.0',
                  'max_min_ratio']

    scatter_matrix(df[attributes], figsize=(20, 16))
    plt.savefig(out, dpi=800)


def _get_train_test(df, test_size=0.1, scale=True):
    '''Splits data into train / test.'''
    x_df = df.select_dtypes(include=['float64'])
    y_df = pd.DataFrame(df['max_min_ratio'], columns=['max_min_ratio'])
    x_df = x_df.drop(['min_mean', 'min_sd', 'max_mean', 'max_sd',
                      'max_min_ratio'], axis=1)

    if scale:
        x_df = _standard_scale(x_df)
        y_df = _standard_scale(y_df)

    return train_test_split(x_df, y_df, test_size=test_size)


def _standard_scale(df):
    '''Perform standard scaling of data.'''
    scaler = StandardScaler()
    return pd.DataFrame(scaler.fit_transform(df), columns=df.columns)


def _rfe(x_train, x_test, y_train, y_test):
    '''Perform linear regression.'''
    lin_reg = LinearRegression()
    rfe = RFE(estimator=lin_reg, n_features_to_select=1, step=1)
    rfe.fit(x_train, y_train)

    for vals in reversed(sorted(zip(rfe.ranking_, x_train.columns))):
        print '\t'.join([str(val) for val in vals])

    _print_result(y_test, rfe.predict(x_test))


def _stab_select(x_train, y_train):
    '''Perform stability selection.'''
    rlasso = RandomizedLasso(alpha=0.025)
    rlasso.fit(x_train, y_train)

    for vals in reversed(sorted(zip(rlasso.scores_, x_train.columns))):
        print '\t'.join([str(val) for val in vals])


def _print_result(y_test, y_pred):
    '''Prints result.'''
    for vals in zip(y_test.values, y_pred):
        print '\t'.join([str(val) for val in vals])

    print 'MSE: ' + str(mean_squared_error(y_test.values, y_pred))


def main(args):
    '''main method.'''
    df = pd.read_csv(args[0])
    results_df = clean(pd.read_csv(args[1]))
    # results_df = results_df.drop(results_df.index[0])
    results_df['max_min_ratio'] = \
        results_df['max_mean'] / results_df['min_mean']

    df = pd.merge(df, results_df, on='variant')
    df.to_csv('out.csv', index=False)

    # print df.describe()
    # df.hist(bins=10, figsize=(20, 15))
    # plt.show()

    # _plot_scatter(df)

    x_train, x_test, y_train, y_test = _get_train_test(df, test_size=0.05)
    _rfe(x_train, x_test, y_train, y_test)
    _stab_select(x_train, y_train)


if __name__ == '__main__':
    main(sys.argv[1:])
