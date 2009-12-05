#       basic_tools.py
#       
#       Copyright 2007 Raja Selvaraj <rajajs@gmail.com>
#       
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

from __future__ import division
import scipy, scipy.signal, scipy.stats.stats
import pylab
import os, commands


"""Some basic signal processing tools"""

def _norm_dot_product(a,b):
    """Will return normalized dot product for two vectors"""
    # TODO: Needs better documentation
    anorm = a/(scipy.dot(a,a)**0.5)
    bnorm = b/(scipy.dot(b,b)**0.5)
    return scipy.dot(anorm,bnorm)

def _ms_to_samples(ms, samplingrate):
    """convert interval in ms to number of samples.
    samplingrate is samples / second"""
    return samplingrate * ms / 1000


def ClassifyBeats(ecg, segmentstart, segmentlength, qtinterval, qrswidth, globalgoodQRS = None, debug=False):
    """
    In the given ECG, find the beats which match an averaged template with less
    than 98% pearson correlation coefficient. For each beat, the segment starts 
    at segmentstart and extends for segmentlength. Returns a vector of same length
    as segmentstart, 1 for good match and 0 for no match. globalgoodQRS is an 
    optional flag vector for good qrs determined 'globally'
    June 17 2008 - changed threshold to 0.90
    """
    
    # Basic checks - do we have a full segment at the end ?
    if len(ecg) < segmentstart[-1] + segmentlength:
        segmentstart = segmentstart[:-1]
    
    # to store debug data
    corr_list = []
    
    # get qt variance
        # and the t wave variance
    t_variance = get_twave_variance(ecg, segmentstart, qtinterval, qrswidth)
    
    # use only 'clean' segments to make template
    if globalgoodQRS != None:
        cleanstart = []
        for segindex in range(len(segmentstart)):
            if globalgoodQRS[segindex] == 1:
                cleanstart.append(segmentstart[segindex])
        cleanstart = scipy.array(cleanstart)
    else:
        cleanstart = segmentstart

    # Make template
    template = scipy.zeros(segmentlength)
    for pt in range(segmentlength):
        template[pt] = scipy.mean(ecg[cleanstart+pt])
    
    # Check each to template and classify 
    localgoodQRS = []
    for start in segmentstart:
        seg = ecg[start:start+segmentlength]
        (corr, pval ) = scipy.stats.stats.pearsonr(seg, template)
        if corr > 0.95:
            localgoodQRS.append(1)
        else:
            localgoodQRS.append(0)
        
        if debug:
            corr_list.append(corr)
            
    ## classify by qt variance
    for beat in range(len(segmentstart)):
        if t_variance[beat] > 3:
            localgoodQRS[beat] = 0
    
    # merge the local and global flags
    goodQRS = scipy.array(localgoodQRS)
    
    if globalgoodQRS != None:
        goodQRS[:] *= globalgoodQRS
    
    if debug:
        return scipy.array(goodQRS), corr_list
    else:
        return scipy.array(goodQRS)

def realign(ecg,QRSonset,QRSwidth, windowSize,SampleRate=1000, debug = False):
    """Syntax newQRSonset = realign(ecg,QRSonset,QRSwidth,windowsize,samplingrate)
    Give ecg and calculated qrsonset
    Realigns by two iterations of maximal normalised dot product
    windowSize - is the initial search window size around the marked QRS to 
    allow for matching - is given in ms
    Raja S"""
    
    #convert windowsize to samples
    windowSize = int(windowSize*SampleRate/1000)
    #Check if first and last beats can be used
    firstbeat = 0
    if QRSonset[0] < windowSize:
        QRSonset = QRSonset[1:] #.pop[0]

    #remove last beat if it is with out a QRSend
    Nbeats = len(QRSonset)
    if QRSonset[-1] + windowSize + QRSwidth > len(ecg):
        lastbeat = Nbeats - 1
    else:
        lastbeat = Nbeats

    #Create first template
    nQRSonset = scipy.zeros(lastbeat,int)  #New QRS onsets
    goodQRSflag = scipy.zeros(lastbeat,int)

    template1 = scipy.zeros(QRSwidth)
    for t in range(QRSwidth):
        template1[t] = scipy.mean(ecg[QRSonset[:lastbeat]+t])
        
    #align and exclude morphologically different beats
    dotproduct = scipy.zeros(2*windowSize + 1)
    dp1 = []  # to get data out for debugging
    dp2 = []
    
    for t in range(lastbeat):
        for point in range(-windowSize,windowSize):
            begin = QRSonset[t]+point
            dotproduct[point+windowSize] = NormDotProduct(template1,ecg[begin:begin+QRSwidth])
        onset = scipy.argmax(dotproduct)
        xx = dotproduct[onset]
        nQRSonset[t] = QRSonset[t] - windowSize - 1 + onset
        if debug:
            dp1.append(xx)
        if xx > 0.6:
            goodQRSflag[t] = 1

    #Create second template using only 'good' beats
    newQRSonset = scipy.zeros(len(nQRSonset),int)
    template2sum = scipy.zeros(QRSwidth)
    template2 = scipy.zeros(QRSwidth)

    for beat in range(lastbeat):
        template2sum[:QRSwidth] = template2sum[:QRSwidth]+(goodQRSflag[beat]*ecg[nQRSonset[beat]:nQRSonset[beat]+QRSwidth])

    template2 = template2sum / sum(goodQRSflag)
    #Align to second template
    #alignWindowwidth = round(6*SampleRate/1000);
    windowSize = int(windowSize/2) #since we have imported division from future
    dotproduct2 = scipy.zeros(2*windowSize + 1)
    
    #print "creating template2"
    
    for t in range(len(nQRSonset)):
        for point in range(-windowSize,windowSize):
            begin = nQRSonset[t]+point
            dotproduct2[point+windowSize+1] = NormDotProduct(template2,ecg[begin:begin+QRSwidth])

        onset = scipy.argmax(dotproduct2)
        xx = dotproduct2[onset]
        #[xx,onset] = max(dotproduct2);
        newQRSonset[t] = nQRSonset[t] - windowSize - 1 + onset
        if debug:
            dp2.append(xx)

        if xx > 0.96:
            goodQRSflag[t] = 1
        else:
            goodQRSflag[t] = 0
            #goodQRSflag[t-1] = 0 # Added this on 20 June 08

    if firstbeat == 2:
        newQRSonset = newQRSonset[2:]
        goodQRSflag = goodQRSflag[2:]
        
    # further modify qrsflag based on CL
    rrinterval = newQRSonset[1:] - newQRSonset[:-1]
    meanrr = int(scipy.mean(rrinterval))
    
    #modify qrsflags by cycle length criterion
    for i in range(len(rrinterval)):
        if rrinterval[i] < 0.85*meanrr:
            #goodQRSflag[i] = 0
            goodQRSflag[i+1] = 0
      

          
    if debug:
        return (newQRSonset,goodQRSflag, dp1, dp2)
    else:
        return (newQRSonset,goodQRSflag)
 

def get_twave_variance(ecg, qrsonset, qtinterval, qrswidth):
    """Quantify variance of T wave amplitude of beat compared to rest.
    Output value if the mean of the absolute z scores of
    t amplitude point by point in the JT interval"""
    
    jtinterval = qtinterval - qrswidth
    Nbeats = len(qrsonset)
    
    ecgmat = scipy.zeros((Nbeats, jtinterval))
    for beat in range(Nbeats):
        ecgmat[beat,:] = ecg[qrsonset[beat]+qrswidth: qrsonset[beat]+qtinterval]
    
    # get t mean and sd at each point on jt interval
    mean_t = []; sd_t = []
    for point in range(jtinterval):
        point_vector = ecgmat[:,point]
        mean_t.append(scipy.mean(point_vector))
        sd_t.append(scipy.std(point_vector))
    
    t_variance = []
    for beat in range(Nbeats):
        zscore  = 0
        for point in range(jtinterval):
            zscore += abs((ecgmat[beat,point] - mean_t[point]) / sd_t[point])
        t_variance.append(zscore / len(qrsonset))    

    return t_variance   
     
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
    
    
def powerSpect(signal):
    """
    Return the power spectrum for the signal
    """
    
    #make sure we use an even number so that the spect contains odd no. of points
    if windowlength%2 == 0:
        Nfft = len(signal)
    else:
        Nfft = len(signal) + 1

    powerspect = scipy.zeros(Nfft/2 + 1)
        
    #remove dc
    signal -= scipy.mean(signal) 
    
    #get the first half of the spectrum    
    spect = scipy.fft(signal,Nfft)[:Nfft/2 + 1]
    
    #get absolute magnitude and scale 
    spect = abs(spect)/Nfft
    
    #powerspect is sq of this
    powerspect = spect**2
    
    #except dc and nyquist, other points have to be multiplied by 2
    powerspect[1:-1] *= 2

    return powerspect

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
        

def vectorAlt(vector):
    """Measure alternans using the spectral
    method in a vector with a sampling rate of 1"""
    timeseries = vector - scipy.mean(vector)  #remove dc

    windowlength = len(timeseries)
    if windowlength % 2 == 0:
        Nfft = windowlength
    else:
        Nfft = windowlength+1


    #get the first half of the spectrum    
    spect = scipy.fft(timeseries,Nfft)[:Nfft/2 + 1]

    #get absolute magnitude and scale it by nbeats
    spect = abs(spect)/Nfft

    #powerspect is sq of this
    powerspect = spect**2

    #except dc and nyquist, other points have to be multiplied by 2
    powerspect[1:-1] *= 2

    #calculate valt and k for point
    altpower = powerspect[-1]
    noise = powerspect[-11:-1]
    meannoise = scipy.mean(noise)
    stdnoise = scipy.std(noise)

    if altpower < meannoise:
        valt = 0
    else:
        valt = scipy.sqrt(altpower - meannoise)

    kvalue = (altpower - meannoise)    / stdnoise

    return valt, kvalue, scipy.sqrt(meannoise), powerspect


#def altMeasure(powerspect):
    #"""
    #Given power spectrum, calculate the k value
    #and Valt
    #"""
    ##alternans power is last value
    #altpower = powerspect[-1]
    
    ##noise is 10 values before 0.5
    #noise = powerspect[-11:-1]
    #meannoise = scipy.mean(noise)
    #stdnoise = scipy.std(noise)
    
    ##valt and k
    #if altpower < meannoise:
        #valt = 0
    #else:
        #valt = scipy.sqrt(altpower - meannoise)
    
    #k = (altpower - meannoise) / stdnoise    
    
    #return valt,k


def breakSegments(total, seglength):
    """
    Return the indices for the start and end
    of window segments when breaking into
    segments with seglength
    """
    segments = []
    startpoints = range(0,total,seglength)
    
    for start in startpoints:
        seg = [start,start+seglength-1]
        if seg[1] > total:
            seg[1] = total
        
        segments.append(seg)
        
    return segments
    

def drawECG(data, savefilename=None):
    """
    Draw a 12 lead ECG with background grid 
    Data is a matrix with 12 columns 
    and atleast 10 seconds of data (at 1000 Hz)
    If savefilename is not given, ecg will be plotted
    """
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
        
    onemm = 40
    onesec = 1000
    onemv = 400
    lenecg = 10*onesec
    htecg = 136*onemm

    #Linethicknesses
    thickbgwidth = 0.4
    thinbgwidth = 0.1
    ecgwidth = 0.8

    pylab.figure()
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
    
if __name__ == "__main__":
    """
    Do all the testing here
    """
    filename = '/media/win_d/Projects_in_progress/finipress/Mcilroy/TWAapace600to500.txt'
    resp = pylab.load(filename,usecols=[12,])
    resp = resp[4000:184000]
    p = countResp(resp,1000)
    
    #pylab.plot(presp)
    #pylab.hold(True)
    #pylab.plot(p[:],presp[p[:]],'xr')
    #pylab.plot(presp)
    #pylab.hold(True)
    #pylab.plot(p[:],presp[p[:]],'xr')
    #pylab.show()    
