import math
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pylab import *

f="data/averagepulse.txt"

averagePulseArray = np.loadtxt(f,delimiter=" ")
t = linspace(0,5,5000)
plotAvg, = plt.plot(t,averagePulseArray,label='Electron Temp')
crystalFig = plt.gcf()
crystalFig.set_size_inches(6,4)
crystalFig.savefig('avgPulse.png',dpi=150)
