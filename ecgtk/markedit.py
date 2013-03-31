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

###################

# Warning: Quickly hacked together as I need it to use for work

#####################
import yaml
import scipy
import wx
import os

from wx.lib.pubsub import Publisher as pub

import matplotlib
from matplotlib.backends.backend_wx import FigureCanvasWx as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure

from io_utils import BardReader
from ecgtk import Cursor
from ecgtk import ECG
##################

ID_LOAD = wx.NewId()
ID_MARK = wx.NewId()
ID_SAVE = wx.NewId()
ID_CHANGE_CHANNEL = wx.NewId()

##################

def _pick_file(msg):
    """
    GUI to pick a single file
    """
    filepicker = wx.FileDialog(None,
                               message=msg, style=wx.OPEN)
    if filepicker.ShowModal() == wx.ID_OK:
        filename = filepicker.GetPath()
        return filename
    else:
        return None

        
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
        self.ECG_READERS = {'cust_numpy' : CustReader,
                            'bard'       : BardReader}

        
    def load_ecg(self, filename = '', format='cust_numpy', marks_filename = ''):
        """
        Should load ECG data and info stored in various formats
        Default format is custom format with numpy and yaml
        Other formats can be specified and are listed in ECG_READERS
        """
        reader = self.ECG_READERS[format](filename)
        self.ecg_data, self.ecg_info = reader.read()
        print 'Loaded ecg'
        if marks_filename != '':
            self.marks = self.load_marks(marks_filename)
            print 'Loaded marks', marks_filename
        else:
            print 'No Marks file'
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
        

    def load_marks(self, filename):
        self.marks = scipy.load(filename)

    def save_marks(self, filename):
        if self.marks == None:
            return
        scipy.save(filename, self.marks)

    def generate_marks(self):
        """
        Generate annotations for the ecg
        """
        print 'starting qrs detection'
        ecg = ECG(self.ecg_data, self.ecg_info)
        qrspeaks = ecg.get_qrspeaks(self.channel1)
        print 'detected %s qrs complexes' %(len(qrspeaks))
        x = len(qrspeaks)
        self.marks = scipy.concatenate(
                                      (qrspeaks.reshape(x, 1),
                                       scipy.ones(x).reshape(x,1)),
                                      1)

        pub.sendMessage("ECG LOADED", (self.ecg_data[:, self.channel1],
                                       self.ecg_data[:, self.channel2],
                                       self.marks))        
        
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

        self.add_menubar()
        
        self.add_toolbar()  
        # http://stackoverflow.com/questions/4740988/add-new-navigate-modes-in-matplotlib
        self.pan_tool  = self.toolbar.FindById(
                         self.toolbar._NTB2_PAN)
        self.zoom_tool = self.toolbar.FindById(
                         self.toolbar._NTB2_ZOOM)

        self.Bind(wx.EVT_CLOSE, self.on_close)


    def on_close(self, event):
        pub.sendMessage("VIEW CLOSED")
        
        
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


    def add_menubar(self):
        """
        Loading and saving will be available in the menu
        """
        menubar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(ID_LOAD, "&Load", "Load ECG and/or annotation")
        file_menu.Append(ID_MARK, "&Annotate", "Generate new annotation")
        file_menu.Append(ID_SAVE, "&Save", "Save annotations")

        channel_menu = wx.Menu()
        channel_menu.Append(ID_CHANGE_CHANNEL, "&Change Chan 1",
                            "Change channel 1")
        
        menubar.Append(file_menu, "&File")
        menubar.Append(channel_menu, "&Channel")
        
        self.SetMenuBar(menubar)
        

    def init_plot(self, ecg1, ecg2, marks):
        """
        Plot will show two ECG leads, each with the marks, if available
        Last row will be rr intervals
        """
        print 'Plotting'
        ecg1, ecg2 = self.remove_overlap(ecg1, ecg2)

        # if zoomed in, maintain x axis zoom
        minx, maxx =  self.axes.get_xlim()

        self.axes.clear()
        self.axes.plot(ecg1, 'b')
        self.axes.plot(ecg2, 'b')
        if marks != None:
            for samp, label in marks:
                self.axes.plot(samp, ecg1[samp], 'ok')
                self.axes.plot(samp, ecg2[samp], 'ok')

        if maxx != 1.0: # which is the default without a plot
            self.axes.set_xlim(minx, maxx)
                
        self.canvas.draw()


    def remove_overlap(self, ecg1, ecg2):
        """
        Adjust baseline of ecg1 so that it doesnt overlap
        with ecg2
        """
        ht_ecg2 = ecg2.max() - ecg2.mean()
        depth_ecg1 = ecg1.mean() - ecg1.min()
        offset = ecg1.mean() - ecg2.mean()
        min_offset = ht_ecg2 + depth_ecg1
        if offset < min_offset:
            ecg1 += min_offset

        return ecg1, ecg2
        
        
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

        self.set_bindings()

    def set_bindings(self):
        cid = self.view.canvas.mpl_connect(
            'button_press_event', self.onclick)
        pub.subscribe(self.init_plot, "ECG LOADED")
        pub.subscribe(self.add_mark_plot, "MARK ADDED")
        pub.subscribe(self.remove_mark_plot, "MARK REMOVED")
        pub.subscribe(self.on_view_closed, 'VIEW CLOSED')

        self.view.Bind(wx.EVT_MENU, self.load_ecg, id=ID_LOAD)
        self.view.Bind(wx.EVT_MENU, self.add_marks, id=ID_MARK)
        self.view.Bind(wx.EVT_MENU, self.save_marks, id=ID_SAVE)
        self.view.Bind(wx.EVT_MENU, self.change_channel, id=ID_CHANGE_CHANNEL)

        
    def load_ecg(self, event):
        ecgfile = _pick_file("Pick the ECG file")
        if ecgfile == None:
            return

        self.ecgfile = ecgfile
        extension = os.path.splitext(ecgfile)[1]
        if extension == '.npy':
            formt = 'cust_numpy'
        elif extension == '.txt':
            formt = 'bard'

        marks_filename = ecgfile.rstrip(extension) + '.ann'
        if not os.path.exists(marks_filename):
            print 'Marks file does not exist ', marks_filename
            marks_filename = ''
            
        self.model.load_ecg(ecgfile, formt, marks_filename)
        

    def add_marks(self, event):
        """
        Generate annotations using Pantomkins for QRS detection
        and add the marks to the view
        """
        self.model.generate_marks()


    def save_marks(self, event):
        """
        Save the current marks to file
        """
        extension = os.path.splitext(self.ecgfile)[1]
        marks_filename = self.ecgfile.rstrip(extension) + '.ann'
        scipy.save(marks_filename, self.model.marks)
        print 'saved to ', marks_filename
        
    def load(self, ecg_data, ecg_info={'samp_freq':1000}, ann=None):
        """
        Load data for the ecg and annotation marks
        and display them
        ecg_data : ndarray, samples x channels
        ecg_info : dict
        ann : n x 2 array, columns are sample no., annot
        """
        self.model.set_data(ecg_data, ecg_info, ann)

    def set_leads(self, lead1, lead2):
        self.model.channel1 = lead1
        self.model.channel2 = lead2


    def change_channel(self, event):
        """
        Allow user to change the channel chosen for first lead
        """
        current_channel1 = self.model.channel1
        try:
            channel_names = self.model.ecg_info['signal_names']
        except KeyError:
            channel_names = self.model.ecg_info['channellabels']

        cc = ChannelChooser(self.view, channel_names, current_channel1)
        if cc.ShowModal() == wx.ID_OK:
            self.model.channel1 = cc.selected_channel
            self.view.init_plot(self.model.ecg_data[:, self.model.channel1],
                                self.model.ecg_data[:, self.model.channel2],
                                self.model.marks)
        else:
            return
        
        
    def init_plot(self, message):
        #lead1, lead2, marks = message.data
        #lead1 += (lead2.max() -  lead2.min())
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

                
    def get_marks(self):
        """
        Return the updated marks
        """
        return self.model.marks


    def on_view_closed(self, message):
        self.view.Destroy()
        return self.model.marks
        
        
class CustReader:
    """
    Reader for the custom ecg format
    """
    def __init__(self, filename):
        self.filename = filename
        self.data = None
        self.info = None
        
    def read(self):
        infofile = os.path.splitext(self.filename)[0] + '.inf'
        data = scipy.load(self.filename)
        info = yaml.load(open(infofile))
        return data, info


class ChannelChooser(wx.Dialog):
    """
    Dialog to choose a channel
    """
    def __init__(self, parent, channel_names, current_channel):
        """
        channel_names is a list of strings
        current_channel is an int
        """
        wx.Dialog.__init__(self, parent, -1, 'Choose new channel')
        self.channel_names = channel_names
        self.current_channel = current_channel

        panel = wx.Panel(self, -1)
        panelsizer = wx.BoxSizer(wx.VERTICAL)
        buttonsizer = wx.BoxSizer(wx.HORIZONTAL)

        # widgets and buttons
        self.choices = wx.RadioBox(panel, -1, "Choose channel",
                                   choices = self.channel_names,
                                   style = wx.RA_VERTICAL)
        self.cancel_button = wx.Button(panel, -1, 'Cancel')
        self.done_button = wx.Button(panel, -1, 'Done')

        self.cancel_button.Bind(wx.EVT_BUTTON, self.cancel)
        self.done_button.Bind(wx.EVT_BUTTON, self.done)

        # add to sizers
        buttonsizer.Add(self.cancel_button, 0, wx.ALL, 10)
        buttonsizer.Add(self.done_button, 0, wx.ALL, 10)
        panelsizer.Add(self.choices, 15, wx.ALL|wx.EXPAND, 2)
        panelsizer.Add(buttonsizer, 1, wx.ALL, 2)

        panel.SetSizer(panelsizer)

        mainsizer = wx.BoxSizer(wx.HORIZONTAL)
        mainsizer.Add(panel, 1, wx.EXPAND, 5)
        mainsizer.Fit(self)
        self.Layout()

        self.SetSize((300, 500))
        

    def cancel(self, event):
        self.EndModal(wx.ID_CANCEL)


    def done(self, event):
        self.selected_channel = self.choices.GetSelection()
        self.EndModal(wx.ID_OK)

        
if __name__ == '__main__':
    from  wfdbtools import rdsamp, rdann
    
    app = wx.App(False)
    controller = Controller(app)

    # testrec = '../samples/format212/100'

    # data, info = rdsamp(testrec, 0, 60)
    # data = data[:, [2,3]] # remove first 2 cols
    
    # ann = rdann(testrec, 'atr', 0, 60)
    # ann = ann[:, [0,2]] # remove second col

    # controller.load(data, info, ann)
    app.MainLoop()

