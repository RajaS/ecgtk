#!/usr/bin/env python

# Raja Selvaraj
from __future__ import division

import wx
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigCanvas
from matplotlib.backends.backend_wxagg import NavigationToolbar2Wx as Toolbar
from matplotlib.widgets import Cursor
from matplotlib.figure import Figure

import numpy

from bard_reader import BardReader


ID_OPEN = wx.NewId()

class EGMPlotter(wx.Frame):
    """Plots electrograms given the data as an array
    of samples x signals"""
    def __init__(self):
        wx.Frame.__init__(self, None, -1, 'Plotter', size=(800,600))
        
        self.plotpanel = PlotPanel(self)
        self.axis = self.plotpanel.fig.gca()

        self._layout()
        self._build_menubar()
        self._set_bindings()

        self.show_cursor()
        self.SHOW_CURSOR = False
        
    def _layout(self):
        mainsizer = wx.BoxSizer()
        mainsizer.Add(self.plotpanel, 1, wx.EXPAND, 5)
        self.SetSizer(mainsizer)
        self.SetSize((800, 600))

    def _build_menubar(self):
        """Build the menu bar"""
        self.menubar = wx.MenuBar()

        file_menu = wx.Menu()
        file_menu.Append(ID_OPEN, "&Open\tCtrl-O","Open file")
        
        self.menubar.Append(file_menu, "&File")
        
        self.SetMenuBar(self.menubar)

        
    def _set_bindings(self):
        """Bindings for menu items and buttons"""
        self.Bind(wx.EVT_MENU, self.select_file, id=ID_OPEN)
        self.plotpanel.canvas.Bind(wx.EVT_KEY_DOWN, self.on_key)
        self.plotpanel.canvas.Bind(wx.EVT_LEFT_DCLICK, self.on_click)


    def on_click(self, event):
        print event.X
        
    def show_cursor(self):
        """display cursor"""
        self.cursor = Cursor(self.axis,  {'linecolor':'red'})
        self.cursor.horizOn = False
        #print dir(self.cursor)

    def toggle_cursor(self):
        """toggle cursor display"""
        print 'caught toggle'
        print self.SHOW_CURSOR
        
        if self.SHOW_CURSOR == False:
            #self.show_cursor()
            self.cursor.visible = True
            self.SHOW_CURSOR = True
            print 'showing'
        else:
            self.cursor.visible = False
            #self.cursor.clear(None)
            self.SHOW_CURSOR = False
            print 'hiding'
        
    def on_key(self, event):
        """key presses for manipulating displayed signals"""
        code = event.GetKeyCode()
        if code == 315: # up arrow
            if self.SELECTED_SIG > 0:
                self.PREV_SELECTION = self.SELECTED_SIG
                self.SELECTED_SIG -= 1
                self.on_selection_change()
        elif code == 317: # down arrow
            if self.SELECTED_SIG < len(self.offsets) - 1:
                self.PREV_SELECTION = self.SELECTED_SIG
                self.SELECTED_SIG += 1
                self.on_selection_change()
        elif code == 61: # +
            self.change_scale(2)
        elif code == 45: # -
            self.change_scale(0.5)
        elif code == 67: # c
            self.toggle_cursor()
        else:
            print code


    def change_scale(self, ratio):
        """Change scale of the selected signal by the ratio"""
        self.scales[self.SELECTED_SIG] *= ratio

        self.xmin, self.xmax = self.axis.get_xlim()
        self.ranges, self.offsets = self.calculate_scaling(self.data)
        self.offsets = self.calculate_offsets(self.scales)
        
        # refresh plots
        self.display_data = self.tru2display(self.data,
                                             self.ranges, self.scales, self.offsets)

        self.plot_data()
        return
            
            
    def on_selection_change(self):
        """change the selected signal.
        Do not change x axis limits"""
        xlim = self.axis.get_xlim()
        self.axis.plot(self.display_data[:, self.PREV_SELECTION], 'black')
        self.axis.plot(self.display_data[:, self.SELECTED_SIG], 'blue')
        self.axis.set_xlim(*xlim)
        
        self.plotpanel.canvas.draw()
        #self.plotpanel.Refresh()
            

    def on_mpl(self, event):
        """catch mouse click in axes"""
        # TODO: work in progress
        if event.inaxes:
            if event.button == 2:
                y = event.ydata

        dist = [offset-y for offset in self.offsets]
        new_selection = dist.index(min(dist))

        event.Skip()
        

    def select_file(self, event):
        """Select file to load"""
        dlg = wx.FileDialog(None, 'Select file to load')
        if dlg.ShowModal() == wx.ID_OK:
            self.datafile = dlg.GetPath()
            self.load_data(self.datafile, 'bard') #TODO: need to be able to specify format
        else:
            return None
        
    def load_data(self, datafile, recorder):
        """datafile is path to file containing data
        format is the recorder used to get the data and 
        defines the reader to use"""
        print 'loading ', datafile
        if recorder == 'bard':
            reader = BardReader(datafile)
        else: # add more later
            pass

        self.data, self.info = reader.data, reader.info

        print 'converting to display data'
        self.init_vals()
        # initial values for xrange
        self.xmin = 0; self.xmax = self.n_samp

        self.ranges, self.offsets = self.calculate_scaling(self.data)

        self.display_data = self.tru2display(self.data,
                                        self.ranges, self.scales, self.offsets)

        self.FIRST_PLOT = True
        
        self.plot_data()

        
    def plot_data(self):
        """Plot the data. Set y limits tight and retain x limits"""
        xlim = self.axis.get_xlim()
        
        self.axis.clear()
        
        self.axis.plot(self.display_data, 'black')
        self.axis.plot(self.display_data[:, self.SELECTED_SIG], 'blue')
        
        self.axis.set_yticks(self.offsets)
        self.axis.set_yticklabels(self.info['signal_names'])

        # redraws maintain x axis limits
        if not self.FIRST_PLOT:
            self.axis.set_xlim(*xlim)
            
        self.axis.set_ylim(self.offsets[-1] - self.scales[-1]/2.0,
                           self.offsets[0] + self.scales[0]/2.0)
        self.plotpanel.canvas.draw()

        self.FIRST_PLOT = False


    def init_vals(self):
        """initialise some values"""
        self.SELECTED_SIG = 0  # one signal is always selected
        self.PREV_SELECTION = 0

        self.n_samp, self.n_sig = self.data.shape
        self.scales = [1] * self.n_sig
        
        
    def calculate_scaling(self, data):
        """calculate current ranges and offsets"""
        ranges = []

        for i in range(self.n_sig):
            sig = data[self.xmin:self.xmax, i]
            ranges.append(max(sig) - min(sig))

        offsets = self.calculate_offsets(self.scales)
        return ranges, offsets


    def calculate_offsets(self, scales):
        """Given the scales of the signals calculate offsets.
        first offset will always be 0.
        subsequent ones must be adjusted so there is no signal overlap"""
        offsets = [0]
        for i in range(len(scales)-1):
            offsets.append(offsets[-1] + (scales[i] + scales[i+1]) / 2.0)

        # invert offsets
        offsets = [offsets[-1] - offset for offset in offsets]
        
        return offsets
    
        
    def tru2display(self, data, ranges, scales, offsets):
        """Convert the array data so that the
        columns can be plotted as signals of equal range
        without overlapping
        """
        display = numpy.copy(data)
        n_samp, n_sig = display.shape

        # scale and offset each signal
        for i in range(n_sig):
            display[:, i] /= (ranges[i] / scales[i])
            display[:, i] += offsets[i]

        return display


        
class PlotPanel(wx.Panel):
    """
    Panel with embedded matplotlib plots
    """
    def __init__(self,parent):
        wx.Panel.__init__(self,parent,-1)
        self.fig = Figure(figsize=(2,2))
        self.canvas = FigCanvas(self,-1,self.fig)

        self.toolbar = Toolbar(self.canvas)
        self.toolbar.Realize()
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.canvas, 1, wx.EXPAND)
        sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        self.SetSizer(sizer)

if __name__ == '__main__':
    app = wx.PySimpleApp()
    plotter = EGMPlotter()
    plotter.Show(True)
    app.MainLoop()
