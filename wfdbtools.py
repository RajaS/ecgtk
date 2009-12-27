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

## Annotation codes
CODEDICT = {
    0 : 'NOTQRS',	# not-QRS (not a getann/putann codedict) */
    1 : 'NORMAL',	# normal beat */
    2 : 'LBBB',	# left bundle branch block beat */
    3 : 'RBBB',	# right bundle branch block beat */
    4 : 'ABERR',	# aberrated atrial premature beat */
    5 : 'PVC',	# premature ventricular contraction */
    6 : 'FUSION',	# fusion of ventricular and normal beat */
    7 : 'NPC',	# nodal (junctional) premature beat */
    8 : 'APC',	# atrial premature contraction */
    9 : 'SVPB',	# premature or ectopic supraventricular beat */
    10 : 'VESC',	# ventricular escape beat */
    11 : 'NESC',	# nodal (junctional) escape beat */
    12 : 'PACE',	# paced beat */
    13 : 'UNKNOWN',	# unclassifiable beat */
    14 : 'NOISE',	# signal quality change */
    16 : 'ARFCT',	# isolated QRS-like artifact */
    18 : 'STCH',	# ST change */
    19 : 'TCH',	# T-wave change */
    20 : 'SYSTOLE',	# systole */
    21 : 'DIASTOLE',	# diastole */
    22 : 'NOTE',	# comment annotation */
    23 : 'MEASURE',	# measurement annotation */
    24 : 'PWAVE',	# P-wave peak */
    25 : 'BBB',	# left or right bundle branch block */
    26 : 'PACESP',	# non-conducted pacer spike */
    27 : 'TWAVE',	# T-wave peak */
    28 : 'RHYTHM',	# rhythm change */
    29 : 'UWAVE',	# U-wave peak */
    30 : 'LEARN',	# learning */
    31 : 'FLWAV',	# ventricular flutter wave */
    32 : 'VFON',	# start of ventricular flutter/fibrillation */
    33 : 'VFOFF',	# end of ventricular flutter/fibrillation */
    34 : 'AESC',	# atrial escape beat */
    35 : 'SVESC',	# supraventricular escape beat */
    36 : 'LINK',	# link to external data (aux contains URL) */
    37 : 'NAPC',	# non-conducted P-wave (blocked APB) */
    38 : 'PFUS',	# fusion of paced and normal beat */
    39 : 'WFON',	# waveform onset */
    #WFON : 'PQ',	# PQ junction (beginning of QRS) */
    40 : 'WFOFF',	# waveform end */
    #WFOFF : 'JPT',	# J point (end of QRS) */
    41 : 'RONT'	# R-on-T premature ventricular contraction */
    }


def rdsamp(record, start=0, end=-1, interval=-1):
    """
    Read signals from a format 212 record from Physionet database.

    Only 2 channel records in format 212 are supported.
    This is the most common record in the
    Physionet database(http://www.physionet.org/physiobank/).
    
    Parameters
    ----------
    record : str
            Full path to record. No extension to be used for record name.
    start  : int, optional
            time to begin in seconds, default 0
    end    : int, optional
            time to end in seconds, defaults to end of record
    interval : int, optional
            interval of data to be read.
            If both interval and end are given, earlier limit is used.

    Returns
    -------
    data : (N, 4) ndarray
          numpy array with 4 columns
          col 1 - Elapsed time in samples
          col 2 - Elapsed time in milliseconds
          col 3,4 - The two signals
          Signal amplitude is in physical units (mV)          
    info : dict
          Dictionary containing header information
          keys :
          'signal_names' - Names of each signal
          'samp_freq' - Sampling freq (samples / second)
          'samp_count' - Total samples in record
          'firstvalues' - First value of each signal
          'gains' - Gain for each signal
          'zerovalues' - Zero value for each signal
    
    """
    # read the header file - output is a dict
    info = rdhdr(record)
    # establish start and end in samples
    start, end = _get_read_limits(start, end, interval, info)
    # read the data
    data = _read_data(record, start, end, info) 
    return data, info

def rdann(record, annotator, start=0, end=-1, types=[]):
    """
    Reads annotations for given record by specified annotator.

    Parameters
    ----------
    record : str
            Full path to record. Record name has no extension.
    annotator : str
            Name of annotator, eg. 'atr'.
            This is the extension for the annotation file.
    start  : int, optional
            time to begin in seconds, default 0
    end    : int, optional
            time to end in seconds, defaults to end of record
    types   : list, optional
            list of annotation types that will be returned.
            Types are input as annotation code (integer from 0 to 49)
            Annotation types not in list will be ignored.
            Default is empty list, which results in all types being read.
            
    Returns
    -------
    data : (N, 3) ndarray
          numpy array with 3 columns
          col 1 - Elapsed time in samples for each annotation.
          col 2 - Elapsed time in seconds for each annotation.
          col 3 - The annotation code.

    """
    # TODO: annotation time is off in the end for 100.atr
    # Error happens after code=61 at annot 1911.
    # last correct annot_time = 546792
    # get header data
    info = rdhdr(record)
    
    annfile = record + '.' + annotator
    fid = open(annfile, 'rb')
    arr = numpy.fromstring(fid.read(), dtype = numpy.uint8).reshape((-1, 2))
    fid.close()

    rows = arr.shape[0]
    annot = []
    annot_time = []
    i = 0

    print rows
    
    while i < rows:
        anntype = arr[i, 1] >> 2
        if anntype == 59:
            annot.append(arr[i+3, 1] >> 2)
            annot_time.append(arr[i+2, 0] + (arr[i+2, 1] << 8) +
                              (arr[i+1, 0] << 16) + arr[i+1, 1] << 24)
            i += 3
        elif anntype in [60, 61, 62, 63]:
            print anntype, i
            print numpy.sum(annot_time)
            hilfe = arr[i, 0] + ((arr[i, 1] & 3) << 8)
            hilfe += hilfe % 2
            i += hilfe / 2
        else:
            annot_time.append(arr[i, 0] + ((arr[i, 1] & 3) << 8))
            annot.append(arr[i, 1] >> 2)
        i += 1
    # annot_time should be total elapsed samples
    annot_time = numpy.cumsum(annot_time)
    annot_time_ms = annot_time / info['samp_freq'] # in seconds
    # limit to requested interval
    start, end = _get_read_limits(start, end, -1, info)
    ann = numpy.array([annot_time, annot_time_ms, annot]).transpose()
    # filter by annot_time in interval
    ann =  ann[start <= ann[:, 0]]
    ann = ann[ann[:, 0] <= end]
    # filter by type
    if types != []:
        ann = ann[numpy.logical_or.reduce([ann[:,2] == x for x in types])]
        #ann = ann[numpy.array([ann[x, 2] in types for x in range(len(ann))])]

    return ann
    
def plot_data(data, info, ann=None):
    """
    Plot the signal with annotations if available.

    Parameters
    ----------
    data : (N, 4) ndarray
         Output array from rdsamp.
    info : dict
         Header information as a dictionary.
         Output from rdsamp
    ann : (N, 2) ndarray, optional
         Output from rdann

    Returns
    -------
    None
    Matplotlib figure is plotted with the signals and annotations.
    
    """
    time = data[:, 1] #in seconds. use data[:, 0] to use sample no.
    sig1 = data[:, 2]
    sig2 = data[:, 3]
    
    pylab.subplot(211)
    pylab.plot(time, sig1, 'k')
    pylab.xticks([])
    pylab.ylabel('%s (mV)' %(info['signal_names'][0]))
    
    pylab.subplot(212)
    pylab.plot(time, data[:, 3], 'k')
    pylab.ylabel('%s (mV)' %(info['signal_names'][1])) 
    pylab.xlabel('Time (seconds)')

    if ann != None:
        # annotation time in samples from start
        ann_x = (ann[:, 0] - data[0,0]).astype('int')
        pylab.plot(ann[:, 1], data[ann_x, 3], 'xr')
        pylab.subplot(211)
        pylab.plot(ann[:, 1], data[ann_x, 2], 'xr')

    pylab.show()

def rdhdr(record):
    """
    Returns the information read from the header file

    Header file for each record has suffix '.hea' and
    contains information about the record and each signal.

    Parameters
    ----------
    record : str
            Full path to record. Record name has no extension.

    Returns
    -------
    info : dict
          Information read from the header as a dictionary.
          keys :
          'signal_names' - Names of each signal
          'samp_freq' - Sampling freq (samples / second)
          'samp_count' - Total samples in record
          'firstvalues' - First value of each signal
          'gains' - Gain for each signal
          'zerovalues' - Zero value for each signal
    
    """
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
    for reading data, determines limits to use.
    info is the dict returned by rdhdr
    end of -1 means end of record.
    If both end and interval are given, choose
    earlier limit of two.
    start and end are returned as samples.
    Example:
    >>> _get_read_limits(0, 2, -1, {'samp_count':100, 'samp_freq':10})
    (0, 20)
    >>> _get_read_limits(0, 2, 3, {'samp_count':100, 'samp_freq':10})
    (0, 20)
    >>> _get_read_limits(0, 4, 2, {'samp_count':100, 'samp_freq':10})
    (0, 20)
    >>> _get_read_limits(0, 105, -1, {'samp_count':100, 'samp_freq':10})
    (0, 100)
    >>> _get_read_limits(-1, -1, -1, {'samp_count':100, 'samp_freq':10})
    (0, 100)
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

    # time columns
    data[:, 0] = numpy.arange(start, end)  # elapsed time in samples
    data[:, 1] = (numpy.arange(samp_to_read) + start) / info['samp_freq'] # in sec
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

def get_annotation_code(code=None):
    """Returns the symbolic definition for the wfdb annotation code.

    See http://www.physionet.org/physiotools/wpg/wpg_31.htm for details.
    Based on ecgcodes.h from wfdb.
    
    Parameters
    ----------
    code : int
           Integer from 0 to 49 (ACMAX).

    Returns
    -------
    Definition : str
                 The definition for the code.

    """
    return CODEDICT[code]

def main():
    """"""
    numpy.set_printoptions(precision=3, suppress=True)
    record  = '/data/Dropbox/programming/ECGtk/samples/format212/100'
    #data, info = rdsamp(record, 10, 20)
    ann = rdann(record, 'atr') #, types=[1])
    #print data
    #print len(data)
    print 'len(ann)', len(ann)
    print ann[1900:1920, :]
    #print info

    #plot_data(data, info, ann)
    
if __name__ == '__main__':
    main()
        
