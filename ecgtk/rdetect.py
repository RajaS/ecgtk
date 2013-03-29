#!/usr/bin/env python

# Design
# ------

# Basic philosophy is to provide a gui for qrs detection / verification and modification in an ecg.
# Should provide for
#   - Input                                              : Read ECG of different formats
#   - QRS detection                                      : Ability to choose of one of different detectors
#                    - Detectors should return qrs marks and labels
#                    - Labels - follow wfdb labels 'http ://www.physionet.org/physiotools/wpg/wpg_31.htm'
#                             - 1 is normal, 5 is PVC
#   - Visualization                                      : See the detected QRS points on a chosen ECG lead
#   - Modification                                       : Delete and add qrs marks. change qrs annotation
#   - Save                                               : Write the marks and labels to file
#                    this will be a numpy array with 2 columns
#                    first col will be sample number of mark
#                    second col will be the label
#   - Load                                               : Load saved marks and labels


#####################
import yaml
import scipy
import wx

from wx.lib.pubsub import Publisher as pub

import matplotlib
from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

from io_utils import BardReader
from ecgtk import Cursor
##################


class Model:
    """
    ECG data - stored as numpy array with as many columns as channels and as many rows as samples
    ECG info - stored as dict
    marks     - stored as numpy array 2 x nmarks, first col is mark, second is label
    """
    def __init__(self):
        self.ecg_data = None
        self.ecg_info = None
        self.marks = None

        self.channel1 = 0 # default channel to show in first position
        self.channel2 = 1 # default channel to show in second pos
        
        #  ECG readers
        ECG_READERS = {'cust_numpy' : CustReader,
                       'bard'       : BardReader}

        
    def load_ecg(self, filename = '', format='cust_numpy'):
        """
        Should load ECG data and info stored in various formats
        Default format is custom format with numpy and yaml
        Other formats can be specified and are listed in ECG_READERS
        """
        reader = ECG_READERS[format]
        self.ecg_data, self.ecg_info = reader.read(filename)
        pub.sendMessage("ECG LOADED", (self.ecg_data[:, self.channel1],
                                       self.ecg_data[:, self.channel2],
                                       self.marks))

    def set_data(self, ecg_data, ecg_info, marks):
        """
        Set the data and info
        """
        self.ecg_data = ecg_data
        self.ecg_info = ecg_info
        self.marks = marks
        pub.sendMessage("ECG LOADED", (self.ecg_data[:, self.channel1],
                                       self.ecg_data[:, self.channel2],
                                       self.marks))

    def add_mark(self, x, label=1):
        """
        Insert new mark into existing marks
        """
        if self.marks == None:
            return

        x = int(x)
        self.marks = scipy.append(self.marks, scipy.array([[x, 1]]), 0)
        self.marks = scipy.sort(self.marks, 0)

        pub.sendMessage("MARK ADDED", (x,
                                       self.ecg_data[x,self.channel1],
                                       self.ecg_data[x,self.channel2]))
        
        
    def remove_mark(self, x):
        """
        Remove existing mark.
        Mark closest to x will be chosen
        """
        print 'removing'
        if self.marks == None:
            return

        closest = scipy.argmin(abs(self.marks[:, 0] - x))
        mark_to_remove = self.marks[closest, 0]
        self.marks = scipy.delete(self.marks, closest, 0)
        pub.sendMessage("MARK REMOVED", (mark_to_remove,
                                       self.ecg_data[mark_to_remove,self.channel1],
                                       self.ecg_data[mark_to_remove,self.channel2]))
        
        
class View(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title='QRS detect')
        self.SetBackgroundColour(wx.NamedColour("WHITE"))
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)

        self.canvas = FigureCanvas(self, -1, self.figure)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()
        
        self.add_toolbar()  
        # http://stackoverflow.com/questions/4740988/add-new-navigate-modes-in-matplotlib
        self.pan_tool  = self.toolbar.FindById(
                         self.toolbar._NTB2_PAN)
        self.zoom_tool = self.toolbar.FindById(
                         self.toolbar._NTB2_ZOOM)
        

        
    def add_toolbar(self):
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.Realize()
        if wx.Platform == '__WXMAC__':
            self.SetToolBar(self.toolbar)
        else:
            # On Windows platform, default window size is incorrect, so set
            # toolbar width to figure width.
            tw, th = self.toolbar.GetSizeTuple()
            fw, fh = self.canvas.GetSizeTuple()
            # By adding toolbar in sizer, we are able to put it at the bottom
            # of the frame - so appearance is closer to GTK version.
            # As noted above, doesn't work for Mac.
            self.toolbar.SetSize(wx.Size(fw, th))
            self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        # update the axes menu on the toolbar
        self.toolbar.update()


    def init_plot(self, ecg1, ecg2, marks):
        """
        Plot will show two ECG leads, each with the marks, if available
        Last row will be rr intervals
        """
        ecg1 += (ecg2.max() -  ecg2.min())

        self.axes.plot(ecg1, 'b')
        self.axes.plot(ecg2, 'b')
        if marks != None:
            for samp, label in marks:
                self.axes.plot(samp, ecg1[samp], 'ok')
                self.axes.plot(samp, ecg2[samp], 'ok')

    def add_mark(self, x, y1, y2):
        """
        New mark to be plotted.
        y1 and y2 are heights for lead1 and lead2
        """
        self.axes.plot(x, y1, 'ok')
        self.axes.plot(x, y2, 'ok')
        self.canvas.draw()
        
    def remove_mark(self, x, y1, y2):
        """
        Remove existing marks from leads 1 and 2
        Draw red cross over to denote remove
        """
        print 'plotting at', x, y1
        self.axes.plot(x, y1, 'xr')
        self.axes.plot(x, y2, 'xr')
        self.canvas.draw()
        

    def _get_xrange(self):
        xmin, xmax = self.axes.get_xbound()
        return xmax - xmin
        
                
    def OnPaint(self, event):
        self.canvas.draw()

        
class Controller:
    def __init__(self, app):
        self.model = Model()
        self.view = View(None)
        self.view.Show()

        cid = self.view.canvas.mpl_connect(
            'button_press_event', self.onclick)
        pub.subscribe(self.init_plot, "ECG LOADED")
        pub.subscribe(self.add_mark_plot, "MARK ADDED")
        pub.subscribe(self.remove_mark_plot, "MARK REMOVED")
        
        
    def load(self, ecg_data, ecg_info={'samp_freq':1000}, ann=None):
        """
        Load data for the ecg and annotation marks
        and display them
        ecg_data : ndarray, samples x channels
        ecg_info : dict
        ann : n x 2 array, columns are sample no., annot
        """
        self.model.set_data(ecg_data, ecg_info, ann)

        
    def init_plot(self, message):
        #lead1, lead2, marks = message.data
        self.view.init_plot(*message.data)

    def add_mark_plot(self, message):
        self.view.add_mark(*message.data)

    def remove_mark_plot(self, message):
        self.view.remove_mark(*message.data)


    def onclick(self, event):
        """
        Click events inside the axis are handled
        only when the pan or zoom buttons are not active
        """
        print 'clicked'
        # button 1 is left, button 3 is right
        # ydata gives the sample number
        if self.view.pan_tool.IsToggled() or \
           self.view.zoom_tool.IsToggled():
            pass # do nothing if pan or zoom are enabled
        else:
            if event.button == 1: # left click
                self.model.add_mark(event.xdata)
            elif event.button == 3: # right click
                self.model.remove_mark(event.xdata)

        
        
class CustReader:
    """
    Reader for the custom ecg format
    """
    def __init__(self):
        self.data = None
        self.info = None
        
    def read(filename):
        infofile = os.path.splitext(filename)[0] + '.inf'
        data = scipy.load(filename)
        info = yaml.load(open(infofile))




if __name__ == '__main__':
    from  wfdbtools import rdsamp, rdann
    
    app = wx.App(False)
    controller = Controller(app)

    testrec = '../samples/format212/100'

    data, info = rdsamp(testrec, 0, 60)
    data = data[:, [2,3]] # remove first 2 cols
    
    ann = rdann(testrec, 'atr', 0, 60)
    ann = ann[:, [0,2]] # remove second col

    controller.load(data, info, ann)
    app.MainLoop()



