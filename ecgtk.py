from __future__ import division
import scipy
import tools

class ECG():
    def __init__(self, ecg, samplingrate = 1000):
        self.ecg = ecg
        self.samplingrate = samplingrate
    

    def _remove_cubic_spline(self, anchors, window):
        """Remove baseline wander by subtracting a cubic spline.
        ecg is vector representing one ecg channel.
        anchors is a vector of isoelectric points (usually qrs onset -20ms)
        anchorwindow is  width of window to use for averaging the amplitude at anchors"""
        windowwidth = tools._ms_to_samples(window, self.samplingrate) / 2
    
        #Do we have enough points before first anchor to use it
        if anchors[0] - window < 0:
            firstbeat = 1
        else:
            firstbeat = 0
    
    #print anchors[:10]
    inputECG -= scipy.mean(inputECG[anchors[:]])  #subtract dc
    #Initialize things
    Nbeats = len(anchors)-firstbeat;
    #x = scipy.zeros(Nbeats)
    y = scipy.zeros(Nbeats)
    #xx = scipy.array(range(int(anchors[firstbeat]),int(anchors[-1])+1))
    xx = scipy.array(range(len(inputECG)))
    x = anchors[firstbeat:]
    
    for i in range(Nbeats):
        y[i] = scipy.mean(inputECG[anchors[firstbeat+i]-anchorwindow:anchors[firstbeat+i]+anchorwindow])
    
    #And do the actual work
    #cleanECG = inputECG;
    tck = scipy.interpolate.splrep(x,y)
    yy = scipy.interpolate.splev(xx,tck)
    
    cleanECG = scipy.zeros(len(inputECG))
    #print "inputecg ",len(inputECG)
    #print "yy ", len(yy)
    #print "anchors(firstbeat) ",anchors[firstbeat]
    
    #for t in range(int(anchors[firstbeat]),int(anchors[-1])):
    for t in range(len(inputECG)):
        #print t
        try:
            cleanECG[t] = inputECG[t] - yy[t]
        except:
            print "failed at t ",t    
    
    return cleanECG    
