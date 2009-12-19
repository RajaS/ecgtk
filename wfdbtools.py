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
import numpy
import pylab

def rdsamp(record, start=0, end=-1, interval=-1, timecol=True):
    """
    Read signals from the specified wfdb record.
    **Only supports format 212 records with 2 signals**
    record - name of record without extension
    start - time to begin in seconds, default 0
    end - time to end in seconds, defaults to end of record
    interval - interval to read.
    If both interval and end are given, earlier limit is used.
    If timecol is True (default), elapsed time in seconds is first column.
    If timecol is False, sample number is given instead of elapsed time.
    """
    # read the header file
    (samp_freq, samp_count, gains,
     zerovalues, firstvalues) = _read_header(record)
    # establish start and end in samples
    start, end = _get_read_limits(start, end, interval, samp_freq, samp_count)
    # read the data
    data = _read_data(record, start, end, samp_freq,
                      zerovalues, firstvalues, gains, timecol)
    return data
    
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
               zerovalues, firstvalues, gains, timecol):
    """Read the binary data for each signal"""
    datfile = record + '.dat'
    samp_to_read = end - start

    # verify against first value in header
    fid = open(datfile, 'rb')
    data = _arr_to_data(numpy.fromstring(fid.read(3),
                        dtype=numpy.uint8).reshape(1,3))
    fid.close()
    if [data[0, 1], data[0, 2]] != firstvalues:
        raise ValueError, 'First value does not match' #TODO: warning
    
    # read into an array with 3 bytes in each row
    fid = open(datfile, 'rb')
    fid.seek(start*3)
    arr = numpy.fromstring(fid.read(3*samp_to_read),
                dtype=numpy.uint8).reshape((samp_to_read, 3))
    fid.close()
    data = _arr_to_data(arr)

    # adjust zerovalue and gain
    data[:, 1] -= zerovalues[0] / gains[0]
    data[:, 2] -= zerovalues[1] / gains[1]
    time = numpy.arange(samp_to_read) *1000 // samp_freq # in ms
    if timecol:
        data[:, 0] = time
    else:
        data[:, 0] = numpy.arange(end)
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
    # data has columns - time, signal1 and signal2
    data = numpy.zeros((arr.shape[0], 3), dtype='int')
    data[:, 1] = (bytes1 << 8) + arr[:, 0] - sign1
    data[:, 2] = (bytes2 << 8) + arr[:, 2] - sign2
    return data

def plot_data(data, ann=None):
    """Plot the signals"""
    time = data[:, 0]
    pylab.subplot(211)
    pylab.plot(time, data[:, 1], 'k')
    pylab.subplot(212)
    pylab.plot(time, data[:, 2], 'k')
    pylab.show()


def test():
    data = rdsamp('/data/Dropbox/programming/ECGtk/samples/format212/100', 1, 7)
    plot_data(data)
    
if __name__ == '__main__':
    test()
        
