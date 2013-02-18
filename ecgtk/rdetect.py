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