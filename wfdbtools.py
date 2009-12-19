#!/usr/bin/env python

"""
Import ECGs in 'format 212' which is the most common format
for ECG data in Physionet.
Only records with 2 channels are supported"""

# Based on rddata.m for matlab written by Robert Tratnig
# available here -
# http://physionet.org/physiotools/matlab/rddata.m

# Author: Raja Selvaraj <rajajs@gmail.com>
# Written in December 2009
# License: GPL

from __future__ import division
import os
import numpy
import pylab

class Record212():
    """WFDB format 212 record"""
    def __init__(self, datafile, samplestart=0, sampleend=0):
        """
        datafile is the path to the dat file
        samplestart and sampleend are in seconds
	"""
        self.datafile = datafile
        recordname = os.path.splitext(datafile)[0]
        self.headerfile = recordname + '.hea'
        self.atrfile = recordname + '.atr'
        self.samplestart = samplestart 
        self.sampleend = sampleend
        # store header information for each signal
        self.gain = []
        self.bitres = []
        self.zerovalue = []
        self.firstvalue = []
        # Load header information
        self.signals = []
        self.load_header()
        # unless samplestart and end are valid,
        # read all the data
        if not self.valid_start_end():
            self.samplestart = 0
            self.sampleend = self.samp_ct

        # read the binary data
        self.read_data()

        # plot data
        self.plot_data()
        
    def valid_start_end(self):
        """are samplestart and end inputs valid"""
        start, end = self.samplestart, self.sampleend
        if start < 0 or end < 0:
            return False
        elif end <= start:
            return False
        else:
            return True
        
    def load_header(self):
        """Load the header data"""
        fi = open(self.headerfile, 'r')
        firstline = fi.readline()
        recordname, signal_ct, samp_rate, samp_ct = firstline.split()

        self.signal_ct = int(signal_ct)
        self.samp_rate = int(samp_rate)
        self.samp_ct = int(samp_ct)

        # only data with 2 signals
        if self.signal_ct != 2:
            raise ValueError, 'Input data must have exactly 2 signals'
        # change samplestart and end to samples
        self.samplestart *= self.samp_rate
        self.sampleend *= self.samp_rate
        
        # load signal attributes
        # TODO: more robust parsing 
        for i in range(self.signal_ct):
            (filename,
             format_name,                 # should be 212
             gain,                        # Integers per mV
             bitres,                      # Bit Resolution
             zerovalue,                   # value of zero point
             firstvalue,                  # First value of signal
             dummy,
             dummy,
             signal_name
             ) = fi.readline().rstrip('\n').split()

            if format_name != '212':
                raise ValueError, 'Not in format 212'
            self.gain.append(int(gain))
            self.bitres.append(int(bitres))
            self.zerovalue.append(int(zerovalue))
            self.firstvalue.append(int(firstvalue))
        fi.close()
            
    def read_data(self):
        """Read the binary data for each signal"""
        fi = open(self.datafile, 'rb')
        # read into an array with 3 bytes in each row
        arr = numpy.fromstring(fi.read(3*self.sampleend),
                    dtype=numpy.uint8).reshape((self.sampleend, 3))
        # bit operations to get the 12-bit data
        second_col = arr[:, 1].astype('int')
        bytes1 = second_col & 15 # bytes belonging to first sample
        bytes2 = second_col >> 4 # belongs to second sample
        sign1 = (second_col & 8) << 9 # sign bit for first sample
        sign2 = (second_col & 128) << 5 # sign bit for second sample
        
        self.data = numpy.zeros((self.sampleend, 2), dtype=numpy.uint16)
        self.data[:, 0] = (bytes1 << 8) + arr[:, 0] - sign1
        self.data[:, 1] += (bytes2 << 8) + arr[:, 2] - sign2
        # verify with first value
        if [self.data[0, 0], self.data[0, 1]] != self.firstvalue:
            raise ValueError, 'First value does not match' #TODO: warning
        # adjust zerovalue and gain
        self.data[:, 0] -= self.zerovalue[0] / self.gain[0]
        self.data[:, 1] -= self.zerovalue[1] / self.gain[1]
        self.time = numpy.arange(self.sampleend) / self.samp_rate


    def plot_data(self):
        """Plot the signals"""
        pylab.subplot(211)
        pylab.plot(self.time, self.data[:, 0], 'k')
        pylab.subplot(212)
        pylab.plot(self.time, self.data[:, 1], 'k')
        pylab.show()



def test():
    record = Record212('100.dat', 0, 10)
    
if __name__ == '__main__':
    test()
        
# case 1
#     M( : , 1) = (M( : , 1) - zerovalue(1));
#     M( : , 2) = (M( : , 2) - zerovalue(1));
#     M = M';
#     M(1) = [];
#     sM = size(M);
#     sM = sM(2)+1;
#     M(sM) = 0;
#     M = M';
#     M = M/gain(1);
#     TIME = (0:2*(SAMPLEEND_2)-1)/sfreq;
# otherwise 
#     disp('Error: Sorting Algorithm For > 2 Signals Not Programmed Yet!');
# end;
# clear A M1H M2H PRR PRL;
# %**************************************************************************
# % Load Attributes Data
# %**************************************************************************
# atrd = fullfile(PATH, ATRFILE);
# fid3 = fopen(atrd,'r');
# A = fread(fid3, [2, inf], 'uint8')';
# fclose(fid3);
# ATRTIME = [];
# ANNOT = [];
# sa = size(A);
# saa = sa(1);
# i = 1;
# while i <= saa
#     annoth = bitshift(A(i,2),-2);
#     if annoth == 59
#         ANNOT = [ANNOT;bitshift(A(i + 3,2),-2)];
#         ATRTIME = [ATRTIME;A(i+2,1) + bitshift(A(i + 2,2),8) +...
#                 bitshift(A(i + 1,1),16) + bitshift(A(i + 1,2),24)];
#         i = i + 3;
#     elseif annoth == 60
#     elseif annoth == 61
#     elseif annoth == 62
#     elseif annoth == 63
#         hilfe = bitshift(bitand(A(i,2),3),8) + A(i,1);
#         hilfe = hilfe + mod(hilfe,2);
#         i = i + hilfe/2;
#     else
#         ATRTIME = [ATRTIME;bitshift(bitand(A(i,2),3),8) + A(i,1)];
#         ANNOT = [ANNOT;bitshift(A(i,2),-2)];
#    end;
#    i = i + 1;
# end;
# ANNOT(length(ANNOT)) = [];                  % Last Line = EOF (= 0)
# ATRTIME(length(ATRTIME)) = [];              % Last Line = EOF
# clear A;
# ATRTIME = (cumsum(ATRTIME))/sfreq;
# ind = find(ATRTIME <= TIME(end));
# ATRTIMED = ATRTIME(ind);
# ANNOT = round(ANNOT);
# ANNOTD = ANNOT(ind);
# %**************************************************************************
# % Manipulate Data So We Only Look At What The User Wants
# %**************************************************************************
# ECG_1_Temp = M(:,1);
# ECG_1 = ECG_1_Temp(SAMPLESTART_2 : SAMPLEEND_2);
# if nosig == 2
#     ECG_2_Temp = M(:,2);
#     ECG_2 = ECG_2_Temp(SAMPLESTART_2 : SAMPLEEND_2);
# end;
# Time_Adjusted = TIME(SAMPLESTART_2 : SAMPLEEND_2);
# %**************************************************************************
# % Display Data
# %**************************************************************************
# figure(1); clf, box on, hold on
# plot(Time_Adjusted, ECG_1,'r');
# if nosig == 2
#     plot(Time_Adjusted, ECG_2,'b');
# end;
# for k = 1:length(ATRTIMED)
#     text(ATRTIMED(k),0,num2str(ANNOTD(k)));
# end;
# xlim([Time_Adjusted(1), Time_Adjusted(end)]);
# xlabel('Time (Seconds)'); ylabel('Voltage (mV)');
# string = ['ECG Signal ',DATAFILE];
# title(string);
# fprintf(1,'\\n$> DISPLAYING DATA FINISHED \n');
# %**************************************************************************
# % Output Data File Into Current Working Directory
# %**************************************************************************
# save(strcat(FILE,'_ECG_',SAMPLESTART,'_',SAMPLEEND) ...
#     , 'ECG_1' , 'ECG_2' , 'Time_Adjusted');
# fprintf(1,'\\n$> ALL FINISHED \n');
# %**************************************************************************
# % End Of Code
# %**************************************************************************
