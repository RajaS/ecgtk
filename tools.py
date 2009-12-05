
from __future__ import division
import scipy, scipy.signal
import pylab
import commands

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
