"""Tests for functions in the wfdbtools module."""

import os
from wfdbtools import rdsamp, rdann, rdhdr

testdir = os.path.abspath('../samples/format212/')
testrecord = os.path.join(testdir, '100')

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

def rdhdr_test():
    """test wfdbtools.rdhdr"""
    for rec in ['100',
                'header_nobells',
                'header_bellsandwhistles']:
        record = os.path.join(testdir, rec)
        info = rdhdr(record)
        assert info['first_values'] == [995.0, 1011.0]
        assert info['gains'] == [200.0, 200.0]
        assert info['signal_names'] == ['MLII', 'V5']
        assert info['units'] == ['mV', 'mV']
        assert info['zero_values'] == [1024.0, 1024.0]

    info = rdhdr(os.path.join(testdir, '7001'))
    assert info['first_values'] == [-53.0, -69.0]
    assert info['gains'] == [100.0, 100.0]
    assert info['samp_count'] == 525000

    # multichannel_header
    record = os.path.abspath('../samples/format16/twa01')
    info = rdhdr(record)
    assert info['signal_count'] == 12
    assert info['signal_names'] ==  ['I', 'II', 'III', 'aVR', 'aVL', 'aVF', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6']
    
    
import nose
nose.main()
