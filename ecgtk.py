from __future__ import division
import scipy
import tools

class ECG():
    def __init__(self, ecg, samplingrate = 1000):
        self.ecg = ecg
        self.samplingrate = samplingrate
    

    def _remove_baseline(self, anchorx, window):
        """Remove baseline wander by subtracting a cubic spline.
        ecg is vector representing one ecg channel.
        anchorx is a vector of isoelectric points (usually qrs onset -20ms)
        window is width of window to use (in ms) for averaging the amplitude at anchors"""
        windowwidth = tools._ms_to_samples(window, self.samplingrate) / 2
    
        #Do we have enough points before first anchor to use it
        if anchorx[0] < windowwidth:
            anchorx = anchorx[1:]

        # subtract dc
        self.ecg -= scipy.mean(self.ecg[anchorx[:]]) 

        # amplitudes for anchors
        anchory = scipy.array([scipy.mean(self.ecg[anchorx[x]-windowwidth:anchorx[x]+windowwidth])
                               for x in range(len(anchorx))])
        # x values for spline that we are going to calculate
        splinex = scipy.array(range(len(self.ecg)))

        # calculate cubic spline fit
        tck = scipy.interpolate.splrep(anchorx, anchory)
        spliney = scipy.interpolate.splev(splinex, tck)

        # subtract the spline
        self.ecg -= spliney
    
        return
    
