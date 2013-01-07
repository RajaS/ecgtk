from __future__ import division
import scipy
import scipy.signal
import pylab
from datetime import datetime
import glob
import os
import matplotlib

import sys
sys.path.append("/data/Dropbox/programming/ECGtk")
import io
#import ecgtk, io

# Running the tests
# run from a terminal "nosetests -v --with-doctest ecgtk.py"
# from io import BardReader

def _norm_dot_product(a,b):
    """Will return normalized dot product for two vectors"""
    # TODO: Needs better documentation
    anorm = a/(scipy.dot(a,a)**0.5)
    bnorm = b/(scipy.dot(b,b)**0.5)
    return scipy.dot(anorm,bnorm)

def _ms_to_samples(ms, samplingrate):
    """convert interval in ms to number of samples.
    samplingrate is samples / second
    >>> _ms_to_samples(500, 1000)
    500
    >>> _ms_to_samples(100, 500)
    50"""
    return int(samplingrate * ms / 1000)

def _samples_to_ms(samples, samplingrate):
    """convert an interval in samples to
    time in ms. samplingrate is samples/second.
    >>> _samples_to_ms(500, 1000)
    500
    >>> _samples_to_ms(50, 500)
    100
    """
    return int(samples * 1000 / samplingrate)

def _format_time_wfdb(ms):
    """convert time in ms to format compatible with rdann.
    This is in the form (hh):mm:ss.sss
    >>> _format_time_wfdb(7322002)
    '02:02:02.002'
    """
    hr, minute = ms//3600000 % 24, ms//60000 % 60
    sec, ms = ms//1000 % 60, ms % 1000
    timeobj = datetime.time(hr, minute, sec, ms*1000) # last val is microsecs
    return timeobj.isoformat()[:-3] # back to ms

def _lfilter_zi(b,a):
    #compute the zi state from the filter parameters. 
    #Based on:
    # Fredrik Gustafsson, Determining the initial states in forward-backward 
    # filtering, IEEE Transactions on Signal Processing, pp. 988--992, April 1996, 
    # Volume 44, Issue 4
    n=max(len(a),len(b))
    zin = (scipy.eye(n-1) - scipy.hstack( (-a[1:n,scipy.newaxis],
                                 scipy.vstack((scipy.eye(n-2), scipy.zeros(n-2))))))
    zid=  b[1:n] - a[1:n]*b[0]
    zi_matrix=scipy.linalg.inv(zin)*(scipy.matrix(zid).transpose())
    zi_return=[]
    #convert the result into a regular array (not a matrix)
    for i in range(len(zi_matrix)):
      zi_return.append(float(zi_matrix[i][0]))

    return scipy.array(zi_return)

def filtfilt(b,a,x):
    """
    Filter with given parameters forward and in reverse to eliminate
    phase shifts.
    In addition, initial state is calculated with lfilter_zi and 
    mirror images of the sample are added at end and beginning to
    remove edge effects.
    Must be a one-dimensional array only.
    """
    #For now only accepting 1d arrays
    ntaps=max(len(a),len(b))
    edge=ntaps*3

    if x.ndim != 1:
        raise ValueError, "Filtfilt is only accepting 1 dimension arrays."

    #x must be bigger than edge
    if x.size < edge:
        raise ValueError, "Input vector needs to be bigger than 3 * max(len(a),len(b)."

    if len(a) < ntaps:
        a=scipy.r_[a,scipy.zeros(len(b)-len(a))]

    if len(b) < ntaps:
        b=scipy.r_[b,scipy.zeros(len(a)-len(b))]

    zi=_lfilter_zi(b,a)

    #Grow the signal to have edges for stabilizing 
    #the filter with inverted replicas of the signal
    s=scipy.r_[2*x[0]-x[edge:1:-1],x,2*x[-1]-x[-1:-edge:-1]]
    #in the case of one go we only need one of the extrems 
    # both are needed for filtfilt

    (y,zf)=scipy.signal.lfilter(b,a,s,-1,zi*s[0])

    (y,zf)=scipy.signal.lfilter(b,a,scipy.flipud(y),-1,zi*y[-1])

    return scipy.flipud(y[edge-1:-edge+1])

def _rms(vector):
    """returns the root mean square
    >>> _rms(scipy.array([1,2,3,4,5]))
    3.3166247903553998
    """
    return scipy.sqrt(scipy.mean(vector**2))

def _zeropad(shortvec, l):
    """Pad the vector shortvec with terminal zeros to length l
    >>> _zeropad(scipy.array([1,2,3,4,5]), 10)
    array([1, 2, 3, 4, 5, 0, 0, 0, 0, 0])
    """
    return scipy.hstack((shortvec, scipy.zeros((l - len(shortvec)), dtype='int')))

def _write_ann(self, qrs_peaks, annfile):
    """Write an annotation file for the QRS onsets in a format
    that is usable with wrann. qrspeaks is in samples"""
    fi = open(annfile, 'w')
    for qrs in qrs_peaks:
        fi.write('%s '*5 + '%s\n' %(_format_time_wfdb(_sample_to_ms(qrs)),
                                    qrs, 'N', 0, 0, 0))
    fi.close()

def get_stim_times(stim, samplingrate):
    """
    stim is a vector of stim recording.
    return a vector of stim times
    """
    blank = 300 * int(samplingrate / 1000)  # min of 300 ms
    threshold = 5000  # 5 mV

    above_threshold = [pt for pt in range(len(stim)) if stim[pt] > threshold]
    stims = [above_threshold[x] for x in range(len(above_threshold))
             if above_threshold[x] - above_threshold[x-1] > blank]

    return stims

def stitch_data(parts_data, parts_info):
    """
    stitch data that comes in multiple parts
    The time in the info is only accurate and overlap is variable
    """
    samplingrate = parts_info[0]['samplingrate']

    # is there overlap
    start_ends = [(info['starttime'], info['endtime']) for
                   info in parts_info]
    start_ends = [(datetime.strptime(start, '%H:%M:%S'),
                   datetime.strptime(end, '%H:%M:%S')) for
                  (start, end) in start_ends]
    for i in range(1, len(start_ends)):
        if not start_ends[i][0] <= start_ends[i-1][1] <= start_ends[i][1]:
            print 'no overlap for recording number %s' %(i+1)
            print start_ends[i-1], start_ends[i]
            #print start_ends
            return None            

    # find overlap
    overlap_rows = []
    for i in range(1, len(start_ends)):
        overlap = (start_ends[i-1][1] - start_ends[i][0]).seconds
        #print 'overlap', overlap

        lastrow = parts_data[i-1][-1, :]
        #print lastrow[:5]
        for row in range(len(parts_data[i])):
            if (parts_data[i][row, :] == lastrow).all() == True:
                #print 'actual overlap', row
                overlap_rows.append(row)


    # stitch data
    datasets = [parts_data[0]]
    for r in range(len(overlap_rows)):
        datasets.append(parts_data[r+1][overlap_rows[r]+1:, :])
    combined_data = scipy.concatenate(datasets)

    # combination info
    combined_info = parts_info[0]
    combined_info['samp_count'] = len(combined_data)
    combined_info['endtime'] = parts_info[-1]['endtime']

    return combined_data, combined_info


def makeMat(ecg,qrsonsets,qrsflags):
    """split ecg to make 2d matrix of the form beats x points.
    bad qrst complexes are replaced by average"""
    
    #convert to millivolts
    #ecg = ecg*1000        # do this later in the calling script if reqd
    
    #Get mean cycle length
    rrintervals = qrsonsets[1:]-qrsonsets[:-1]
    meanrr = int(scipy.mean(rrintervals))
        
    # do we have one rr after last qrsonset?
    # am not changing qrsflags, but last qrsflag will become unpaired
    if qrsonsets[-1] + meanrr > len(ecg):
        qrsonsets = qrsonsets[:-1]
    
    #segment into a 2d matrix
    Nbeats = len(qrsonsets)
    ecgmat = scipy.zeros((Nbeats,meanrr))
    
    #for l in range(leads):
    for i in range(Nbeats):
        if qrsflags[i] == 1:
            ecgmat[i,:] = ecg[qrsonsets[i]:qrsonsets[i]+meanrr]
    
    #get mean qrst and insert where qrsflag is 0
    
    #for l in range(leads):
    meanqrst = scipy.zeros(meanrr)
    sumqrst = scipy.zeros(meanrr)
    
    for i in range(meanrr):
        sumqrst[i] += scipy.sum(ecgmat[:,i])
    meanqrst = sumqrst/scipy.sum(qrsflags)
      
    # if there is atleast one good qrs, replace bad by avg good
    if sum(qrsflags) > 0:
        for i in range(Nbeats):
            if qrsflags[i] == 0:
                #print meanqrst.shape
                ecgmat[i,:] = meanqrst[:]
            
    return (ecgmat,meanqrst)


def altMeasure(ecgmat,windowbegin,windowlength,qrswidth,qtinterval,meanrr):
    """Measure alternans for the ecglead defined in the defined window 
    given the matrix of qrst points created by makemat. Calculate for 
    whole cycle, return values for each point and also overall value
    in the qt interval alone"""
    
    #calculate power spectrum for each column
    #make sure we use an even number so that the spect contains odd no. of points
    if windowlength % 2 == 0:
        Nfft = windowlength
    else:
        Nfft = windowlength+1

    beats, meanrr = ecgmat.shape  #mean rr is the second dim of ecgmat

    powerspect = scipy.zeros((Nfft/2 + 1, meanrr))
    kvector = scipy.zeros(meanrr)
    valtvector = scipy.zeros(meanrr)
    
    for i in range(meanrr):
        
        timeseries = ecgmat[windowbegin:windowbegin+windowlength,i]
        timeseries -= scipy.mean(timeseries)  #remove dc

        #get the first half of the spectrum    
        spect = scipy.fft(timeseries,Nfft)[:Nfft/2 + 1]
        
        #get absolute magnitude and scale it by nbeats
        spect = abs(spect)/Nfft
        
        #powerspect is sq of this
        powerspect[:,i] = spect**2
        
        #except dc and nyquist, other points have to be multiplied by 2
        powerspect[1:-1,i] *= 2
        
        #calculate valt and k for point
        altpower = powerspect[-1,i]
        noise = powerspect[-11:-1,i]
        meannoise = scipy.mean(noise)
        stdnoise = scipy.std(noise)
        
        if altpower < meannoise:
            valtvector[i] = 0
        else:
            valtvector[i] = scipy.sqrt(altpower - meannoise)

        kvector[i] = (altpower - meannoise)    / stdnoise
        
    #calculate aggregate power spectrum for st interval only
    avgpowerspect = scipy.zeros(Nfft/2+1)
    
    for i in range(int(Nfft/2)+1):
        avgpowerspect[i] = scipy.mean(powerspect[i,qrswidth:qtinterval])
        
    altpower = avgpowerspect[-1]
    noise = avgpowerspect[-11:-1]
    meannoise = scipy.mean(noise)
    stdnoise = scipy.std(noise)        
    
    valt = scipy.sqrt(altpower - meannoise).real  # only the real part
    k = (altpower - meannoise)    / stdnoise
        
    return (k, valt, meannoise, kvector, valtvector, avgpowerspect)


def analyseTWA(k,valt):
    """Provide automated analysis of TWA results. Input is the matrix
    of k values and the matrix of Valt"""
    
    segs,leads = k.shape
    posleads = []
    maxvalt = 0
    maxlead = -1
    
    for l in range(leads):
        for s in range(segs):
            if k[s,l] >= 3:
                if valt[s,l] >= 1.9:
                    if l not in posleads:
                        posleads.append(l)
                    if valt[s,l] > maxvalt:
                        maxvalt = valt[s,l]
                        maxlead = l
    
    return (posleads,maxvalt,maxlead) 


class QRSDetector():
    """
    """
    def __init__(self, ecgdata, samplingrate=1000):
        """
        - 'ecgdata' : array - points x leads in case of multiple leads
                      or a vector in case of a single
        - 'samplingrate' : samples per second, default 1000
        """
        try:
            self.data = scipy.array(ecgdata, dtype='float')
        except ValueError, msg:
            raise ValueError("Invalid format for ecg data - %s" %(msg))
            
        self.samplingrate = samplingrate

        # convert vector to column array
        if len(self.data.shape) == 1:
            self.data = scipy.array([self.data]).transpose()

        self.points, self.leads = self.data.shape
        if len(self.data.shape) > 1 and self.leads > self.points:
            raise ValueError("ECG data has more columns than rows")

        # we need atleast 8 seconds of data (for initializing buffers)
        if self.points < self.samplingrate * 8:
            raise ValueError("Length of ECG is less than 8 seconds")
        
    def qrs_detect(self, qrslead=0):
         """Detect QRS onsets using modified PT algorithm
         """
         # If ecg is a vector, it will be used for qrs detection.
         # If it is a matrix, use qrslead (default 0)
         if len(self.data.shape) == 1:
             self.raw_ecg = self.data
         else:
             self.raw_ecg = self.data[:,qrslead]

         # butterworth bandpass filter 5 - 15 Hz
         self.filtered_ecg = self._bpfilter(self.raw_ecg)
         # differentiate
         self.diff_ecg  = scipy.diff(self.filtered_ecg)
         # take absolute value (was square in original PT implementation)
         self.abs_ecg = abs(self.diff_ecg)
         # integrate 
         self.int_ecg = self._mw_integrate(self.abs_ecg)
         
         # Construct buffers with last 8 values 
         self._initializeBuffers(self.int_ecg)

         # collect all unique local peaks in the integrated ecg
         peaks = self.peakDetect(self.int_ecg)

         # classify each peak as QRS or noise
         self.checkPeaks(peaks, self.int_ecg)


         # compensate for delay during integration
         self.QRSpeaks -= 40 * (self.samplingrate / 1000)
         
         return self.QRSpeaks

    def qrs_detect_multiple_leads(self, leads=[]):
        """Use multiple leads for qrs detection.
        Leads to use may be given as list of lead indices.
        Default is to use all leads"""
        # leads not specified, switch to all leads
        if leads == []:
            leads = range(self.leads)

        # qrs detection for each lead
        qrspeaks = []
        for lead in leads:
            qrspeaks.append(self.qrs_detect(lead))

        # DEBUG
        print "length of qrs in different channels"
        print [len(x) for x in qrspeaks]

        # zero pad detections to match lengths
        maxlength = max([len(qrspeak_lead) for qrspeak_lead in
                         qrspeaks])
        for lead in range(len(qrspeaks)):
            qrspeaks[lead] = self._zeropad(qrspeaks[lead], maxlength)

        #DEBUG
        print "max length ", maxlength
        print [len(x) for x in qrspeaks]
        
        qrspeaks_array = scipy.array(qrspeaks).transpose()
        self.QRSpeaks = self.multilead_peak_match(qrspeaks_array)
        return self.QRSpeaks

    def _zeropad(self, shortvec, l):
        """Pad the vector shortvec with terminal zeros to length l"""
        return scipy.hstack((shortvec, scipy.zeros((l - len(shortvec)), dtype='int')))
        
    def write_ann(self, annfile):
        """Write an annotation file for the QRS onsets in a format
        that is usable with wrann"""
        fi = open(annfile, 'w')
        for qrspeak in self.QRSpeaks:
            fi.write('%s %s %s %s %s %s\n' %(self._sample_to_time(qrspeak), qrspeak, 'N', 0, 0, 0))
        fi.close()

    def _sample_to_time(self, sample):
        """convert from sample number to a string representing
        time in a format required for the annotation file.
        This is in the form (hh):mm:ss.sss"""
        time_ms = int(sample*1000 / self.samplingrate)
        hr, min, sec, ms = time_ms//3600000 % 24, time_ms//60000 % 60, \
                           time_ms//1000 % 60, time_ms % 1000
        timeobj = datetime.time(hr, min, sec, ms*1000) # last val is microsecs
        return timeobj.isoformat()[:-3] # back to ms
         
    def visualize_qrs_detection(self, savefilename = False):
        """Plot the ecg at various steps of processing for qrs detection.
        Will not plot more than 10 seconds of data.
        If filename is input, image will be saved"""
        ecglength = len(self.raw_ecg)
        ten_seconds = 10 * self.samplingrate
        
        if ecglength > ten_seconds:
            segmentend = ten_seconds
        elif ecglength < ten_seconds:
            segmentend = ecglength

        segmentQRSpeaks = [peak for peak in self.QRSpeaks if peak < segmentend]

        pylab.figure()
        pylab.subplot(611)
        pylab.plot(self.raw_ecg[:segmentend])
        pylab.ylabel('Raw ECG', rotation='horizontal')
        pylab.subplot(612)
        pylab.plot(self.filtered_ecg[:segmentend])
        pylab.ylabel('Filtered ECG',rotation='horizontal')
        pylab.subplot(613)
        pylab.plot(self.diff_ecg[:segmentend])
        pylab.ylabel('Differential',rotation='horizontal')
        pylab.subplot(614)
        pylab.plot(self.abs_ecg[:segmentend])
        pylab.ylabel('Squared differential',rotation='horizontal')
        pylab.subplot(615)
        pylab.plot(self.int_ecg[:segmentend])
        pylab.ylabel('Integrated', rotation='horizontal')
        pylab.subplot(616)
        pylab.hold(True)
        pylab.plot(self.raw_ecg[:segmentend])
        pylab.plot(segmentQRSpeaks, self.raw_ecg[segmentQRSpeaks], 'xr')
        pylab.hold(False)
        pylab.ylabel('QRS peaks', rotation='horizontal')

        if savefilename:
            pylab.savefig(savefilename)
        else:
            pylab.show()
        
    def _initializeBuffers(self, ecg):
        """Initialize the 8 beats buffers using values
        from the first 8 one second intervals        
        """
        onesec = self.samplingrate
        # signal peaks are peaks in the 8 segments
        self.signal_peak_buffer = [max(ecg[start*onesec:(start+1)*onesec])
                                                  for start in range(8)]
        self.noise_peak_buffer = [0] * 8
        self.rr_buffer = [1] * 8
        self._updateThreshold()
        
    def _updateThreshold(self):
        """Calculate threshold based on amplitudes of last
        8 signal and noise peaks"""
        noise = scipy.mean(self.noise_peak_buffer)
        signal = scipy.mean(self.signal_peak_buffer)
        self.threshold = noise + 0.3125 * (signal - noise)

    def peakDetect(self, ecg):
        """Determine local maxima that are larger than others in
        adjacent 200
        """
        # list all local maxima
        peak_indices = [i for i in range(1,len(ecg)-1)
                     if ecg[i-1] < ecg[i] > ecg[i+1]]
        peak_amplitudes = [ecg[peak] for peak in peak_indices]

        # restrict to peaks that are larger than anything else 200 ms
        # on either side
        unique_peaks = []
        minimumRR = self.samplingrate * 0.2

        # start with first peak
        peak_candidate_index = peak_indices[0]
        peak_candidate_amplitude = peak_amplitudes[0]

        # test successively against other peaks
        for peak_index, peak_amplitude in zip(peak_indices, peak_amplitudes):
            # if new peak is less than minimumRR away and is larger,
            # it becomes candidate
            if peak_index - peak_candidate_index <= minimumRR and\
                                  peak_amplitude > peak_candidate_amplitude:
                peak_candidate_index = peak_index
                peak_candidate_amplitude = peak_amplitude

            # if new peak is more than 200 ms away, candidate is promoted to
            # a unique peak and new peak becomes candidate
            elif peak_index - peak_candidate_index > minimumRR:
                unique_peaks.append(peak_candidate_index)
                peak_candidate_index = peak_index
                peak_candidate_amplitude = peak_amplitude

            else:
                pass

        return unique_peaks

    def checkPeaks(self, peaks, ecg):
        """Check the given peaks one by one according to
        thresholds that are constantly updated"""
        #amplitudes = [ecg[peak] for peak in peaks]
        self.QRSpeaks = [0] # will remove zero later
        
        # augment the peak list with the last point of the ecg
        peaks += [len(ecg)-1]
        
        for index in range(len(peaks)-1):
            peak = peaks[index]
            amplitude = ecg[peak]
            amp_ratio = amplitude / self.threshold
            # accept as QRS if larger than threshold
            # slope in raw signal +-30% of previous slopes - not implemented
            if amp_ratio > 1:
                self.acceptasQRS(peak, amplitude)

            # reject if less than half threshold
            elif amp_ratio < 0.5:
                self.acceptasNoise(peak, amplitude)
                
            # acccept as qrs if higher than half threshold,
            # but is 360 ms after last qrs and
            # next peak is more than 1.5 rr intervals away
            # just abandon it if there is no peak before or after
            else:
                meanrr = scipy.mean(self.rr_buffer)
                lastQRS_to_this_peak = (peak - self.QRSpeaks[-1]) / self.samplingrate
                lastQRS_to_next_peak = peaks[index+1] - self.QRSpeaks[-1]

                if lastQRS_to_this_peak > 0.36 and lastQRS_to_next_peak > 1.5 * meanrr:
                    self.acceptasQRS(peak, amplitude)
                else:
                    self.acceptasNoise(peak, amplitude)

        self.QRSpeaks = scipy.array(self.QRSpeaks[1:])
        return

    def acceptasQRS(self, peak, amplitude):
        self.QRSpeaks.append(peak)

        self.signal_peak_buffer.pop(0)
        self.signal_peak_buffer.append(amplitude)

        if len(self.QRSpeaks) > 1:
            self.rr_buffer.pop(0)
            self.rr_buffer.append(self.QRSpeaks[-1] - self.QRSpeaks[-2])

    def acceptasNoise(self, peak, amplitude):
        self.noise_peak_buffer.pop(0)
        self.noise_peak_buffer.append(amplitude)
            
    def _mw_integrate(self, ecg):
        """
        Integrate the ECG signal over a defined
        time period. 
        """
        # window of 80 ms - better than using a wider window
        window_length = int(80 * (self.samplingrate / 1000))
        int_ecg = scipy.zeros_like(ecg)
        cs = ecg.cumsum()
        int_ecg[window_length:] = (cs[window_length:] -
                                   cs[:-window_length]) / window_length
        int_ecg[:window_length] = cs[:window_length] / scipy.arange(
                                                   1, window_length + 1)
        return int_ecg

    def _bpfilter(self, ecg):
         """Bandpass filter the ECG with a bandpass setting of
         5 to 15 Hz"""
         # relatively basic implementation for now
         Nyq = self.samplingrate / 2
         wn = [5/ Nyq, 15 / Nyq]
         b,a = scipy.signal.butter(2, wn, btype = 'bandpass')
         # TODO: filtfilt should be implemented here
         return filtfilt(b,a,ecg)

    def multilead_peak_match(self, peaks):
        """Reconcile QRS detections from multiple leads.
        peaks is a matrix of peak_times x leads.
        If the number of rows is different,
        pad shorter series with zeros at end"""
        ms90 = 90 * self.samplingrate / 1000
        Npeaks, Nleads = peaks.shape
        current_peak = 0
        unique_peaks = []

        while current_peak < len(peaks):
            all_values = peaks[current_peak, :]
            outer = all_values.max()
            outerlead = all_values.argmax()
            inner = all_values.min()
            innerlead = all_values.argmin()

            #
            near_inner = sum(all_values < inner + ms90)
            near_outer = sum(all_values > outer - ms90)

            #all are within 90 ms
            if near_inner == near_outer == Nleads:
                unique_peaks.append(int(scipy.median(all_values)))
                current_peak += 1

            # max is wrong
            elif near_inner > near_outer:
                peaks[current_peak+1:Npeaks, outerlead] = peaks[current_peak:Npeaks-1, outerlead]
                peaks[current_peak, outerlead] = scipy.median(all_values)
                # do not change current peak now

            # min is wrong
            elif near_inner <= near_outer:
                peaks[current_peak:Npeaks-1, innerlead] = peaks[current_peak+1:Npeaks, innerlead]
                peaks[-1, innerlead] = 0

        return unique_peaks


class Cursor:
    """
    Cursor that can be used for measurements or
    marking points on a pylab plot
    """
    def __init__(self, ax, cursor='vertical'):
        # cursor can be vertical, horizontal or cross
        self.ax = ax
        self.lx = None; self.ly = None
        self.cursor = cursor

        if cursor != 'vertical':
            self.lx = ax.axhline(color='k')  # the horiz line
        if cursor != 'horizontal':
            self.ly = ax.axvline(color='k')  # the vert line

        pylab.connect('motion_notify_event', self.mouse_move)
        pylab.connect('button_press_event', self.mouse_click)

    def mouse_move(self, event):
        if not event.inaxes: return
        x, y = event.xdata, event.ydata

        # update the line positions
        if self.lx:
            self.lx.set_ydata(y )
        if self.ly:
            self.ly.set_xdata(x )

        pylab.draw()

    def mouse_click(self, event):
        if not event.inaxes: return
        x, y = event.xdata, event.ydata

        if self.cursor == 'vertical':
            self.ax.axvline(x=x, color='r')
        elif self.cursor == 'horizontal':
            self.ax.axhline(y=y, color='r')
        #self.txt.set_text('clicked at %1.2f, %1.2f' %(x,y))


class ECG():
    def __init__(self, data, info = {'samplingrate': 1000}):
        """
        data is a numpy matrix, either single column or multiple
        info is a dict
        units should be in mv
        """
        self.data = data
        self.samplingrate = info['samplingrate']
        self.qrsonsets = None

    def remove_baseline(self, anchorx, window, lead=0):
        """
        Remove baseline wander by subtracting a cubic spline.
        anchorx is a vector of isoelectric points (usually qrs onset -20ms)
        window is width of window to use (in ms) for averaging the amplitude at anchors
        """
        ecg = self.data[:, lead]                    
        windowwidth = _ms_to_samples(window, self.samplingrate) / 2
        #Do we have enough points before first anchor to use it
        if anchorx[0] < windowwidth:
            anchorx = anchorx[1:]
        # subtract dc
        ecg -= scipy.mean(ecg[anchorx[:]]) 
        # amplitudes for anchors
        # window is zero, no averaging
        if windowwidth == 0:
            anchory = scipy.array([ecg[x] for x in anchorx])
        # or average around the anchor
        else:
            anchory = scipy.array([scipy.mean(ecg[x-windowwidth:x+windowwidth])
                      for x in anchorx])
        # x values for spline that we are going to calculate
        splinex = scipy.array(range(len(ecg)))
        # calculate cubic spline fit
        tck = scipy.interpolate.splrep(anchorx, anchory)
        spliney = scipy.interpolate.splev(splinex, tck)
        # subtract the spline
        ecg -= spliney

        self.data[:, lead] = ecg

        return ecg


    def realign(self, qrsonset, qrswidth, windowsize, samplerate, lead):
        """
        Give ecg and calculated qrsonset
        Realigns by two iterations of maximal normalised dot product
        windowSize - is the initial search window size around the marked QRS to 
        allow for matching - is given in ms
        Raja S
        """
        ecg = self.data[:, lead]      
        #convert windowsize to samples
        windowsize = int(windowsize*samplerate/1000)
        #Check if first and last beats can be used
        firstbeat = 0
        if qrsonset[0] < windowsize:
            qrsonset = qrsonset[1:] #.pop[0]

        #remove last beat if it is with out a QRSend
        nbeats = len(qrsonset)
        if qrsonset[-1] + windowsize + qrswidth > len(ecg):
            lastbeat = nbeats - 1
        else:
            lastbeat = nbeats

        #Create first template
        nQRSonset = scipy.zeros(lastbeat,int)  #New QRS onsets
        goodQRSflag = scipy.zeros(lastbeat,int)

        template1 = scipy.zeros(qrswidth)
        for t in range(qrswidth):
            template1[t] = scipy.mean(ecg[qrsonset[:lastbeat]+t])
            
        #align and exclude morphologically different beats
        dotproduct = scipy.zeros(2*windowsize + 1)
        
        for t in range(lastbeat):
            for point in range(-windowsize,windowsize):
                begin = qrsonset[t]+point
                dotproduct[point+windowsize] = _norm_dot_product(template1,
                                                ecg[begin:begin+qrswidth])
            onset = scipy.argmax(dotproduct)
            xx = dotproduct[onset]
            nQRSonset[t] = qrsonset[t] - windowsize - 1 + onset
            if xx > 0.6:
                goodQRSflag[t] = 1

        #Create second template using only 'good' beats
        newQRSonset = scipy.zeros(len(nQRSonset),int)
        template2sum = scipy.zeros(qrswidth)
        template2 = scipy.zeros(qrswidth)

        for beat in range(lastbeat):
            template2sum[:qrswidth] = template2sum[:qrswidth]+ \
                            (goodQRSflag[beat]*\
                            ecg[nQRSonset[beat]:nQRSonset[beat]+qrswidth])

        template2 = template2sum / sum(goodQRSflag)
        #Align to second template
        #alignWindowwidth = round(6*SampleRate/1000);
        windowsize = int(windowsize/2) 
        dotproduct2 = scipy.zeros(2*windowsize + 1)
        
        for t in range(len(nQRSonset)):
            for point in range(-windowsize,windowsize):
                begin = nQRSonset[t]+point
                dotproduct2[point+windowsize+1] = _norm_dot_product(template2,
                                            ecg[begin:begin+qrswidth])

            onset = scipy.argmax(dotproduct2)
            xx = dotproduct2[onset]
            newQRSonset[t] = nQRSonset[t] - windowsize - 1 + onset

            if xx > 0.96:
                goodQRSflag[t] = 1
            else:
                goodQRSflag[t] = 0

        if firstbeat == 2:
            newQRSonset = newQRSonset[2:]
            goodQRSflag = goodQRSflag[2:]
            
        # further modify qrsflag based on CL
        rrinterval = newQRSonset[1:] - newQRSonset[:-1]
        meanrr = int(scipy.mean(rrinterval))
        
        #modify qrsflags by cycle length criterion
        for i in range(len(rrinterval)):
            if rrinterval[i] < 0.85*meanrr:
                goodQRSflag[i+1] = 0

        return newQRSonset, goodQRSflag


    def get_qrspeaks(self, qrslead):
        """
        Using pan tomkins method detect qrs onsets
        currently only qrs peak is detected
        """
        det = QRSDetector(self.data[:, qrslead], self.samplingrate)
        qrsonsets = det.qrs_detect()
        return qrsonsets

    def get_wavelimits(self, qrspeaks, leads=range(12)):
        """
        Given qrspeaks / point on qrs,
        interactively, obtain qrs onset, end and tend
        leads is a list of the indices of ECG leads
        """
        ax = pylab.subplot(111)
        ax.set_title("Pick QRS onset, end and T end")
        #ax = matplotlib.pyplot.axes()
        meanrr = int(scipy.mean(qrspeaks[1:] - qrspeaks[:-1]))
        onems = int(self.samplingrate / 1000)
        r = qrspeaks[int(len(qrspeaks) * 2/3)]  # choose a beat 2/3 of way
        
        start = r - 200 * onems  # 400 ms before
        end = start + meanrr

        for l in leads:
            ax.plot(self.data[start:end, l])

        cursor = Cursor(ax)

        pts = pylab.ginput(3)
        q, s, t = [pt[0] for pt in pts]
        #pylab.show()
        qrsonsets = qrspeaks + int(q - 200 * onems)
        qrsends = qrspeaks + int(s - 200 * onems)
        tends = qrspeaks + int(t - 200 * onems)
        return qrsonsets, qrsends, tends


    def write_ann(self, annfile):
        """Write an annotation file for the QRS onsets in a format
        that is usable with wrann"""
        fi = open(annfile, 'w')
        for qrspeak in self.QRSpeaks:
            fi.write('%s '*4 + '%s\n' %(self._sample_to_time(qrspeak), qrspeak, 'N', 0, 0, 0))
        fi.close()


    def drawECG(self, start=0, leads=range(12), savefilename=None):
        """
        Draw a 12 lead ECG with background grid 
        start is time of recording to start from in seconds
        first 12 leads are used by default, else leads can be specified
        If savefilename is not given, ecg will be plotted
        """
        if self.data.shape[1] < 12:
            raise ValueError, 'Less than 12 leads available'
        if self.data.shape[0] / self.samplingrate < 10:
            raise ValueError, 'Less than 10 seconds of data available'

        data = self.data[:, leads]

        ################################################################################
        #
        #    Draw the background    
        #
        ################################################################################

        #time scale  - 1  ms/px  (1mm = 40 px at 25 mm/s)
        #amp scale -  2.5 mcV/px (1mm = 40 px at 100 mcv/mm) 
        #image width = 10 seconds = 10000 px
        #image height = 4 strips
        #strip height = 34 mm (3.4 mV)
        
        onemm = 40 * int(self.samplingrate / 1000)
        onesec = self.samplingrate
        onemv = (400 * self.samplingrate) / (1000 * 1000) # correct for microV
        lenecg = 10*onesec
        htecg = 136*onemm

        #Linethicknesses
        thickbgwidth = 0.4
        thinbgwidth = 0.1
        ecgwidth = 0.8

        ecgfig = pylab.figure()
        #thick horizontal lines
        for horiz in range(0,htecg,5*onemm):
            pylab.plot([0,lenecg],[horiz,horiz],'r',linewidth=thickbgwidth)

        #thick vertical lines
        for vert in range(0,lenecg,5*onemm):
            pylab.plot([vert,vert],[0,htecg],'r',linewidth=thickbgwidth)

        #thin horizontal lines
        for horiz in range(0,htecg,onemm):
            pylab.plot([0,lenecg],[horiz,horiz],'r',linewidth=thinbgwidth)

        #thin vertical lines
        for vert in range(0,lenecg,onemm):
            pylab.plot([vert,vert],[0,htecg],'r',linewidth=thinbgwidth)
            
        ################################################################################
        #
        #    Draw the ECG    
        #
        ################################################################################
        startplot = 0
        stripcenter = 17 #in mm
        striplength = int(62.5*onemm) # in px (2.5 seconds)
        rhythmlead = 1

        horizcenters = (((scipy.array([0,-1,-2,-3]))*2*stripcenter) - 17) *onemm

        #prepare data
        for lead in range(12):
            #center horizontally
            data[:,lead] -= scipy.mean(data[:,lead])
        #rescale    
        data *= onemv    

        #column 1
        for lead in range(3):
            pylab.plot(range(striplength),\
                       data[startplot:striplength+startplot,lead]-horizcenters[3-lead],\
                       'k',linewidth = ecgwidth)

        #column 2
        for lead in range(3,6):
            pylab.plot(range(striplength,2*striplength),\
                       data[striplength+startplot:striplength*2+startplot,lead]-horizcenters[6-lead],\
                       'k',linewidth = ecgwidth)
            
        #column 3
        for lead in range(6,9):
            pylab.plot(range(2*striplength,3*striplength),\
                       data[striplength*2+startplot:striplength*3+startplot,lead]-horizcenters[9-lead],\
                       'k',linewidth = ecgwidth)

        #column 4
        for lead in range(9,12):
            pylab.plot(range(3*striplength,4*striplength),\
                       data[striplength*3+startplot:striplength*4+startplot,lead]-horizcenters[12-lead],\
                       'k',linewidth = ecgwidth)

        #rhythm strip
        pylab.plot(range(4*striplength),\
                   data[startplot:4*striplength,rhythmlead]-horizcenters[0],\
                   'k',linewidth = ecgwidth)

        ################################################################################
        #
        #    Labels
        #
        ################################################################################
        labels = ['I','II','III','aVR','aVL','aVF','V1','V2','V3','V4','V5','V6']
        xoffset = 20
        yoffset = -250
        labelx = [xoffset]*3               +    [xoffset+striplength]*3 +\
                 [xoffset+2*striplength]*3 +    [xoffset+3*striplength]*3
        labely = list(yoffset - horizcenters[3:0:-1])*4

        for labelct in range(12):
            pylab.text(labelx[labelct],labely[labelct],labels[labelct],fontsize=8)

        pylab.text(labelx[0],yoffset-horizcenters[0],labels[rhythmlead],fontsize=8)

        #pylab.axis('off')
        pylab.setp(pylab.gca(),xticklabels=[])
        pylab.setp(pylab.gca(),yticklabels=[])
        pylab.axis([0,lenecg,0,htecg])
        
        if not savefilename:
            pylab.show()
        else:
            pylab.savefig(savefilename, dpi=300)

            #if possible, crop with imagemagick
            try:
                commands.getoutput("mogrify -trim '%s'" %(savefilename))
            except:
                pass

        # clear the figure
        # Otherwise subsequent plots overlap
        ecgfig.clf()
                

def test_remove_baseline():
    """test for remove_baseline function
    """
    testsignal = scipy.sin(scipy.arange(0,2*scipy.pi,0.01))

    npoints = len(testsignal)
    anchors = range(0,len(testsignal), len(testsignal)//8)
    window = 0

    ecg = ECG(testsignal)
    rms_with_baseline = _rms(ecg.data)
    ecg.remove_baseline(anchors, window)
    rms_without_baseline = _rms(ecg.data)
    assert rms_without_baseline / rms_with_baseline < 0.01


def test():
    from io import BardReader
    f = '/data/Dropbox/work/jipmer_research/post_MI_risk/patient_data/first_case/avpace90_3.txt'
    #f = '/data/Dropbox/work/jipmer_research/post_MI_risk/patient_data/first_case/nsr.txt'    
    br = BardReader(f)
    data, info = br.read()

    print 'loaded data', data.shape

    print info

    ecg = ECG(data, info)
    qrspeaks = ecg.get_qrspeaks(7)

    print 'found qrspeaks', len(qrspeaks)
 
    stim = get_stim_times(ecg.data[:, 13], 2000)

    # ecg.remove_baseline(ecg.qrsonsets-240, 20)

    # pylab.plot(ecg.data[:,7], 'r')
    # for q in ecg.qrsonsets:
    #     pylab.plot(q-240, 400, 'xr')
    # pylab.show()
    # for r in ecg.qrsonsets:
    #     pylab.plot(r, 400, 'xr')
    qrsonsets, qrsends, tends = ecg.get_wavelimits(qrspeaks)

    pylab.plot(ecg.data[:,7])
    for r in qrsonsets:
        pylab.plot(r, 10, 'xr')
    for s in qrsends:
        pylab.plot(s, 10, 'ok')
    for t in tends:
        pylab.plot(t, 10, 'or')
    pylab.show()


def stitch_test():
    from io import BardReader
    filename = "/data/tmp/twa/avpace90_1.txt"
    # find all files in the set
    dirname, basefilename = os.path.split(filename)
    parts = glob.glob(filename[:-5] + '*')
    parts.sort()
    print parts

    # get times
    parts_data = []
    parts_info= []    
    for p in parts:
        br = io.BardReader(p)
        data, info = br.read()
        parts_data.append(data)
        parts_info.append(info)

    stitch_data(parts_data, parts_info)


if __name__ == '__main__':
    stitch_test()
