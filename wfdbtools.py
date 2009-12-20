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
    Read signals from the specified wfdb record.
    **Only supports format 212 records with 2 signals**
    record - name of record without extension
    start - time to begin in seconds, default 0
    end - time to end in seconds, defaults to end of record
    interval - interval to read.
    If both interval and end are given, earlier limit is used.
    Elapsed time in samples is in first column of output.
    Elapsed time in milliseconds is second column.
    """
    # TODO: implement ability to read 
    # read the header file
    (samp_freq, samp_count, gains,
     zerovalues, firstvalues) = _read_header(record)
    # establish start and end in samples
    start, end = _get_read_limits(start, end, interval, samp_freq, samp_count)
    # read the data
    data = _read_data(record, start, end, samp_freq,
                      zerovalues, firstvalues, gains)
    return data


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
    (samp_freq, samp_count, gains,
                   zerovalues, firstvalues) = _read_header(record)
    
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
    annot_time_ms = annot_time * 1000 // samp_freq
    # limit to requested interval
    start, end = _get_read_limits(start, end, -1, samp_freq, samp_count)
    ann = numpy.array([annot_time, annot_time_ms, annot]).transpose()
    # filter by annot_time in interval
    ann =  ann[start <= ann[:, 0]]
    ann = ann[ann[:, 0] <= end]
    return ann
    
def plot_data(data, ann=None):
    """Plot the signals"""
    time = data[:, 1] # used data[:, 2] to use sample no.
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
    gains = []; zerovalues = []; firstvalues = []
    fid = open(headerfile, 'r')
    firstline = fid.readline()
    recordname, signal_count, samp_freq, samp_count = firstline.split()
    # Number of signal must be exactly 2
    if int(signal_count) != 2:
        raise ValueError, 'Input data must have exactly 2 signals'
    
    # load signal attributes
    # TODO: more robust parsing 
    for i in range(int(signal_count)):
        (filename,
         format_name,                 # should be 212
         gain,                        # Integers per mV
         bitres,                      # Bit Resolution
         zerovalue,                   # value of zero point
         firstvalue,                  # First value of signal
         dummy,
         dummy,
         signal_name
         ) = fid.readline().rstrip('\n').split()

        if format_name != '212':
            raise ValueError, 'Not in format 212'
        
        gains.append(int(gain))
        zerovalues.append(int(zerovalue))
        firstvalues.append(int(firstvalue))
    fid.close()
    return (int(samp_freq), int(samp_count), gains,
            zerovalues, firstvalues)

def _get_read_limits(start, end, interval, samp_freq, samp_count):
    """
    Given start time, end time and interval
    for reading data, determine limits to use.
    samp_count is number of samples in record.
    end of -1 means end of record.
    If both end and interval are given, choose
    earlier limit of two.
    """
    start *= samp_freq
    end *= samp_freq
    
    if start < 0:         # If start is negative, start at 0
        start = 0
    if end < 0:           # if end is negative, use end of record
        end = samp_count
    if end < start:       # if end is before start, swap them
        start, end = end, start
    interval_end = start + interval * samp_freq # end determined by interval
    if interval_end < start:
        interval_end = samp_count
    end = min(end, interval_end, samp_count) # use earlier end
    return start, end
            
def _read_data(record, start, end, samp_freq,
               zerovalues, firstvalues, gains):
    """Read the binary data for each signal"""
    datfile = record + '.dat'
    samp_to_read = end - start

    # verify against first value in header
    fid = open(datfile, 'rb')
    data = _arr_to_data(numpy.fromstring(fid.read(3),
                        dtype=numpy.uint8).reshape(1,3))
    fid.close()
    if [data[0, 2], data[0, 3]] != firstvalues:
        warnings.warn('First value from dat file does not match value in header')
    
    # read into an array with 3 bytes in each row
    fid = open(datfile, 'rb')
    fid.seek(start*3)
    arr = numpy.fromstring(fid.read(3*samp_to_read),
                dtype=numpy.uint8).reshape((samp_to_read, 3))
    fid.close()
    data = _arr_to_data(arr)

    # adjust zerovalue and gain
    data[:, 2] -= zerovalues[0] / gains[0]
    data[:, 3] -= zerovalues[1] / gains[1]
    time = numpy.arange(samp_to_read) *1000 // samp_freq # in ms
    data[:, 0] = numpy.arange(start, end) + start
    data[:, 1] = time + (start * 1000 // samp_freq)
    return data

def _arr_to_data(arr):
    """From the numpy array read from file
    using bit level operations,
    extract the 12-bit data"""
    second_col = arr[:, 1].astype('int')
    bytes1 = second_col & 15 # bytes belonging to first sample
    bytes2 = second_col >> 4 # belongs to second sample
    sign1 = (second_col & 8) << 9 # sign bit for first sample
    sign2 = (second_col & 128) << 5 # sign bit for second sample
    # data has columns - samples, time(ms), signal1 and signal2
    data = numpy.zeros((arr.shape[0], 4), dtype='int')
    data[:, 2] = (bytes1 << 8) + arr[:, 0] - sign1
    data[:, 3] = (bytes2 << 8) + arr[:, 2] - sign2
    return data

def test():
    data = rdsamp('/data/Dropbox/programming/ECGtk/samples/format212/100', 1, 10)
    ann = rdann('/data/Dropbox/programming/ECGtk/samples/format212/100', 'atr', 1, 10)

    plot_data(data, ann)
    
    
if __name__ == '__main__':
    test()
        
