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
import matplotlib
from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

from io import BardReader
##################


class Model:
    """
    ECG data - stored as numpy array with as many columns as channels and as many rows as samples
    ECG info - stored as dict
    marks    - stored as numpy array 2 x nmarks, first col is mark, second is label
    """
    def __init__(self):
        ecg_data = None
        ecg_info = None
        marks = None

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
        data, info = reader.read(filename)
        return data, info



class View(wx.Frame):
    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title='QRS detect')
        self.SetBackgroundColour(wx.NamedColour("WHITE"))
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)

        #self.axes.plot(t, s)
        self.canvas = FigureCanvas(self, -1, self.figure)

        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.SetSizer(self.sizer)
        self.Fit()

        self.add_toolbar()  

        
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


    def plot(ecg1, ecg2, marks=None):
        """
        Plot will show two ECG leads, each with the marks, if available
        Last row will be rr intervals
        """
        if marks == None:
            self.axes.plot(ecg1)

        
    def OnPaint(self, event):
        self.canvas.draw()
        
class Controller:
    def __init__(self, app):
        self.model = Model()
        self.view = View(None)

        self.view.Show()

        
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
    app = wx.App(False)
    controller = Controller(app)
    app.MainLoop()



