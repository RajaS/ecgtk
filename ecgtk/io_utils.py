#!/usr/bin/env python

# Raja Selvaraj

# License: GPL

from __future__ import division
import numpy


class BardReader():
    """Read data and header information from files exported
    from Bard EP system"""
    def __init__(self, datafile):
        """datafile is the full path to the exported file"""
        self.datafile = datafile

    def read(self):
        # read header
        with open(self.datafile) as fi:
            self.header = self.get_header(fi)

        # extract some header information
        info, amp_range = self.parse_header_info(self.header)
        
        with open(self.datafile) as fi:
            data = self.read_data(fi)

        # convert data values to microV
        data = self.in_microV(data, info, amp_range)
        info['units'] = 'microV'

        return data, info


    def rows(self, nrows=0):
        """
        Generator returning the data as rows
        Returns upto nrows rows or all data
        """
        with open(self.datafile) as fi:
            self.header = self.get_header(fi)
        info, amp_range = self.parse_header_info(self.header)

        with open(self.datafile) as fi:
            # dump header
            for i in range(len(self.header)):
                h = fi.readline()
                
            if nrows == 0:
                nrows = info['samp_count']
            for i in xrange(nrows):
                yield fi.readline()
            

    def in_microV(self, data, info, amp_range):
        """
        convert the data values in microV
        amp_range is a list, 
            each value is the analog range in V for each channel
        """
        # 2 - extend range on either side
        # 1000 - convert to microV
        # 16 - bit depth
        for chan in range(info['signal_count']):
            m = 2 * amp_range[chan] * 1000 / 2 ** 16
            data[:, chan] = data[:, chan] * m
        return data


    def get_header(self, fi):
        """fi is file object to exported text file.
        Return the header information"""
        header = []
        for l in fi:
            if not l.startswith('[Data]'):
                header.append(l)

            else:
                header.append('Last line')
                return header


    def parse_header_info(self, header):
        """Extract required information from the header
           units - unit of amplitude
           start_time - time of recording start
           end_time - time of recording end
          'signal_count' - Number of signals
          'signal_names' - Names of each signal
          'samp_freq' - Sampling freq (samples / second)
          'samp_count' - Total samples in record
        """
        info = {}
        info['signal_names'] = []

        # channel specific info
        amp_range = []  # ampl range

        for line in header:
            if line.startswith('Channels exported'):
                info['signal_count'] = int(line.split(':')[1].rstrip('\r\n'))
            elif line.startswith('Samples per channel'):
                info['samp_count'] = int(line.split(':')[1].rstrip('\r\n'))
            elif line.startswith('Start time'):
                info['start_time'] = line.lstrip('Start time:').rstrip('\r\n')
            elif line.startswith('End time'):
                info['end_time'] = line.lstrip('End time:').rstrip('\r\n')
            elif line.startswith('Sample Rate'):
                info['samp_freq'] = int(line.split(':')[1].rstrip('Hz\r\n'))

            # extract channel labels
            elif line.startswith('Label'):
                info['signal_names'].append(line.split(':')[1].strip())
            elif line.startswith('Range'): 
                amp_range.append(float(line.split(':')[1].rstrip('mv \r\n')))

            else:
                continue

        return info, amp_range


    def read_data(self, fi):
        """Extract data into numpy array"""
        data = numpy.loadtxt(fi, dtype='float', delimiter=',',
                             skiprows=len(self.header))
        
        
        return data

    
def test():
    f =  '../samples/bard/nsr.txt'
    br = BardReader(f)

    data, info = br.read()
    
    print 'info'
    print info

    print data.shape
    print data

    assert data.shape == (info['samp_count'], info['signal_count'])
    
if __name__ == '__main__':
    test()
