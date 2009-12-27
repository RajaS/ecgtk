"""Tests for functions in the wfdbtools module."""

import os
from wfdbtools import rdsamp, rdann

testrecord = os.path.abspath('samples/format212/100')


def rdsamp_test():
    """test wfdbtools.rdsamp"""
    data, info = rdsamp(testrecord)
    assert info['samp_freq'] == 360
    assert info['signal_names'] == ['MLII', 'V5']
    assert len(data) == 650000
    assert data[0, 2:].tolist() == [-0.145, -0.065]
    assert data[-1, 2:].tolist() == [-1.28, 0.]
    assert int(data[-1, 1]) == 1805
    
    data, info = rdsamp(testrecord, 10, 20)
    assert len(data) == 3600
    assert data[0,2:].tolist() == [-0.39, -0.275]
    assert data[-1, 2:].tolist() == [-0.42, -0.4]

    data, info = rdsamp(testrecord, 10, interval=10)
    assert len(data) == 3600
    assert data[0,2:].tolist() == [-0.39, -0.275]
    assert data[-1, 2:].tolist() == [-0.42, -0.4]

def rdann_test():
    """test wfdbtools.rdann"""
    ann = rdann(testrecord, 'atr')
    assert len(ann) == 2274
    assert ann[0, :].tolist() == [18., 0.05, 28 ]
    assert ann[-1, 0] == 649991.
    
import nose
nose.main()
