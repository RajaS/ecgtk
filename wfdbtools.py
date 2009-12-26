#!/usr/bin/env python

"""
Provide some of the wfdb IO tools in pure python.
I will try to keep names and usage similar to the
original wfdb applications for simplicity.
"""
# Author: Raja Selvaraj <rajajs@gmail.com>
# Written in December 2009
# License: GPL

# rdsamp based on rddata.m for matlab written by Robert Tratnig
# available at http://physionet.org/physiotools/matlab/rddata.m

from __future__ import division
import warnings
import numpy
import pylab

def rdsamp(record, start=0, end=-1, interval=-1):
    """
    Read signals from a format 212 record from Physionet database.

    Only 2 channel records in format 212 are supported.
    This is the most common record in the
    Physionet database(http://www.physionet.org/physiobank/).
    
    Parameters
    ----------
    record : str
            name of record without extension
    start  : int, optional
            time to begin in seconds, default 0
    end    : int, optional
            time to end in seconds, defaults to end of record
    interval : int, optional
            interval of data to be read.
            If both interval and end are given, earlier limit is used.

    Returns
    -------
    data : (N, 3) ndarray
          numpy array with 4 columns
          col 1 - Elapsed time in samples
          col 2 - Elapsed time in milliseconds
          col 3,4 - The two signals
          Signal amplitude is in physical units (mV)          
    info : dict
          Dictionary containing header information

    Example
    -------
    >> data, info = rdsamp('samples/format212/100', 0, 10)
    
    """
    # read the header file - output is a dict
    info = _read_header(record)
    # establish start and end in samples
    start, end = _get_read_limits(start, end, interval, info) 
    # read the data
    data = _read_data(record, start, end, info) 
    return data, info

def rdann(record, annotator, start=0, end=-1):
    """Read the annotation for given record by the annotator.
    record is path to record, no extension.
    annotator is the name of annotator, eg. 'atr'.
    Optional start, end and interval are in seconds
    and specify interval when annot should be read.
    If both end and interval are given, the earlier
    of the two limits is chosen.
    Array that is returned has time in samples in first column,
    time in ms in second column
    and annotations code in third column
    """
    # get header data
    info = _read_header(record)
    
    annfile = record + '.' + annotator
    fid = open(annfile, 'rb')
    arr = numpy.fromstring(fid.read(), dtype = numpy.uint8).reshape((-1, 2))
    fid.close()

    rows = arr.shape[0]
    annot = []
    annot_time = []
    i = 0
    while i < rows:
        anntype = arr[i, 1] >> 2
        if anntype == 59:
            annot.append(arr[i+3, 1] >> 2)
            annot_time.append(arr[i+2, 0] + (arr[i+2, 1] << 8) +
                              (arr[i+1, 0] << 16) + arr[i+1, 1] << 24)
            i += 3
        elif anntype in [60, 61, 62, 63]:
            hilfe = arr[i, 0] + ((arr[i, 1] & 3) << 8)
            hilfe += hilfe % 2
            i += hilfe / 2
        else:
            annot_time.append(arr[i, 0] + ((arr[i, 1] & 3) << 8))
            annot.append(arr[i, 1] >> 2)
        i += 1
    # annot_time should be total elapsed samples
    annot_time = numpy.cumsum(annot_time)
    annot_time_ms = annot_time * 1000 // info['samp_freq']
    # limit to requested interval
    start, end = _get_read_limits(start, end, -1, info)
    ann = numpy.array([annot_time, annot_time_ms, annot]).transpose()
    # filter by annot_time in interval
    ann =  ann[start <= ann[:, 0]]
    ann = ann[ann[:, 0] <= end]
    return ann
    
def plot_data(data, info, ann=None):
    """Plot the signals"""
    # TODO: check if ann has been given
    time = data[:, 1] # use data[:, 2] to use sample no.
    pylab.subplot(211)
    pylab.plot(time, data[:, 2], 'k')
    pylab.plot(ann[:, 1], data[ann[:, 2], 2], 'xr')
    pylab.subplot(212)
    pylab.plot(time, data[:, 3], 'k')
    pylab.plot(ann[:, 1], data[ann[:, 2], 3], 'xr')
    pylab.show()

def _read_header(record):
    """Read the headerfile for the record"""
    headerfile = record + '.hea'
    info = {}
    info['gains'] = []; info['zerovalues'] = []
    info['firstvalues'] = []; info['signal_names'] = []
    fid = open(headerfile, 'r')
    firstline = fid.readline()
    recordname, signal_count, samp_freq, samp_count = firstline.split()
    info['samp_freq'] = int(samp_freq)
    info['samp_count'] = int(samp_count)
    # Number of signal must be exactly 2
    if int(signal_count) != 2:
        raise ValueError, 'Input data must have exactly 2 signals'
    
    # load signal attributes
    # TODO: more robust parsing 
    for line_count in range(int(signal_count)):
        (dummy,                       # filename
         format_name,                 # should be 212
         gain,                        # Integers per mV
         dummy,                       # Bit Resolution
         zerovalue,                   # value of zero point
         firstvalue,                  # First value of signal
         dummy,
         dummy,
         signal_name
         ) = fid.readline().rstrip('\n').split()

        if format_name != '212':
            raise ValueError, 'Not in format 212'
        
        info['gains'].append(int(gain))
        info['zerovalues'].append(int(zerovalue))
        info['firstvalues'].append(int(firstvalue))
        info['signal_names'].append(signal_name)
    fid.close()
    return info

def _get_read_limits(start, end, interval, info):
    """
    Given start time, end time and interval
    for reading data, determine limits to use.
    samp_count is number of samples in record.
    end of -1 means end of record.
    If both end and interval are given, choose
    earlier limit of two.
    """
    start *= info['samp_freq']
    end *= info['samp_freq']
    
    if start < 0:         # If start is negative, start at 0
        start = 0
    if end < 0:           # if end is negative, use end of record
        end = info['samp_count']
    if end < start:       # if end is before start, swap them
        start, end = end, start
    interval_end = start + interval * info['samp_freq'] # end det by interval
    if interval_end < start:
        interval_end = info['samp_count']
    end = min(end, interval_end, info['samp_count']) # use earlier end
    return start, end
            
def _read_data(record, start, end, info):
    """Read the binary data for each signal"""
    datfile = record + '.dat'
    samp_to_read = end - start

    # verify against first value in header
    fid = open(datfile, 'rb')
    data = _arr_to_data(numpy.fromstring(fid.read(3),
                        dtype=numpy.uint8).reshape(1,3))
    fid.close()

    if [data[0, 2], data[0, 3]] != info['firstvalues']:
        warnings.warn(
            'First value from dat file does not match value in header')
    
    # read into an array with 3 bytes in each row
    fid = open(datfile, 'rb')
    fid.seek(start*3)
    arr = numpy.fromstring(fid.read(3*samp_to_read),
                dtype=numpy.uint8).reshape((samp_to_read, 3))
    fid.close()
    data = _arr_to_data(arr)

    # adjust zerovalue and gain
    data[:, 2] = (data[:, 2] - info['zerovalues'][0]) / info['gains'][0]
    data[:, 3] = (data[:, 3] - info['zerovalues'][1]) / info['gains'][1]
    time = numpy.arange(samp_to_read) *1000 // info['samp_freq'] # in ms
    data[:, 0] = numpy.arange(start, end) + start
    data[:, 1] = time + (start * 1000 // info['samp_freq'])
    return data

def _arr_to_data(arr):
    """From the numpy array read from the dat file
    using bit level operations, extract the 12-bit data"""
    second_col = arr[:, 1].astype('int')
    bytes1 = second_col & 15 # bytes belonging to first sample
    bytes2 = second_col >> 4 # belongs to second sample
    sign1 = (second_col & 8) << 9 # sign bit for first sample
    sign2 = (second_col & 128) << 5 # sign bit for second sample
    # data has columns - samples, time(ms), signal1 and signal2
    data = numpy.zeros((arr.shape[0], 4), dtype='float')
    data[:, 2] = (bytes1 << 8) + arr[:, 0] - sign1
    data[:, 3] = (bytes2 << 8) + arr[:, 2] - sign2
    return data

def test():
    """Run some tests"""
    record  = '/data/Dropbox/programming/ECGtk/samples/format212/100'
    data, info = rdsamp(record, 0, 10)
    ann = rdann(record, 'atr', 0, 10)
    print data
    print ann
    print info

    plot_data(data, info, ann)
    
    
if __name__ == '__main__':
    test()
        
