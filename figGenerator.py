import matplotlib.pylab as plt 
import sys
import pickle as pk
import numpy as np
import math
import analyzeGotran as ao
import subprocess as sb
from subprocess import PIPE
import shlex
import time
import ScriptRunner as SR
import matplotlib.pyplot as plt; plt.rcdefaults()
import datetime
import numpy as np
import scipy as sp
import scipy.fftpack
import pandas as pd
from spectrum import *
from scipy.interpolate import spline

## Experimental data 
noP2X7 = np.array([100.00,104.00,168.80,128.53,121.60,118.93,123.20,116.27,108.80,109.87,107.73,110.67,
                   105.33,109.07,102.93,101.33,105.60,104.00,106.13,104.00,113.87,105.33,104.53,103.20,
                   109.33,106.67,103.20,102.67,104.27,102.40]) 
time1 =  np.array([0.00,0.42,0.61,0.98,1.10,1.33,1.44,1.86,2.05,2.08,2.35,2.54,2.73,2.92,3.33,3.64,3.79,3.98,4.09,4.17,
                   4.58,4.89,5.15,5.42,5.68,5.95,6.14,6.52,6.63,6.86])

lowATP = np.array([100.00,100.53,136.80,178.93,174.93,130.13,107.47,100.53,100.27,100.80,100.27,
                   101.60,101.87,100.53,100.80,100.80,101.33,100.53,102.93])
time2  = np.array([0.00,0.38,0.49,1.02,1.21,1.67,2.16,2.77,3.11,4.24,4.32,4.66,4.89,5.08,5.30,5.49,5.91,6.10,6.40])

All    = np.array([100.00,100.00,175.73,166.13,156.27,161.60,153.60,162.93,155.47,164.53,162.40,163.47,154.40,160.27,
                   151.73,161.60,136.27,157.33,137.07,131.73,134.13,134.13,157.87,131.73,146.93,128.80,129.60,129.07,
                   131.47,151.73,125.87,133.33,125.87,149.60,126.40])
time3  = np.array([0.00,0.42,0.61,0.72,1.06,1.17,1.33,1.48,1.93,2.01,
                   2.16,2.20,2.42,2.50,2.69,2.80,3.18,3.41,3.60,3.83,
                   3.94,4.09,4.28,4.70,5.27,5.80,5.91,6.10,6.29,6.52,
                   7.08,7.50,7.92,8.26,8.79])

## WT 
data1 = {}
data1   = SR.gotranMicroglia(sim_time      = 600,
                            ATP           = 1000,
                            output_name   = 'test1',
                            ode_file_name = 'p2xp2yMigration32',
                            rhop2x4       = 1,
                            rhop2x7       = 1,
                            rhop2yc       = 1,
                            data_name2    = 'CaER',
                            data_name3    = 'Cap2y',
                            EleSwitch     = 0,
                            DegSwitch     = 1,
                            removePickle  = 1)