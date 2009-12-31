ECGtk
=====
Provide a toolkit with some simple functions for processing electrocardiograms (ECGs).

wfdbtools.py
============
Provide some functions to read wfdb ECG files from Physionet.

    >> from wfdbtools import rdsamp, rdann, plot_data

    # record is a format 212 record frm physiobank
    >> record  = '/samples/format212/100'

    # Read in the data from 0 to 10 seconds
    >> data, info = rdsamp(record, 0, 10)

    # and the annotation
    >> ann = rdann(record, 'atr', 0, 10)

    # Plot the data and the mark the annotations
    >> plot_data(data, info, ann)

Detailed documentation at http://rajas.github.com/ecgtk/


