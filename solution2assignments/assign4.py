#!/usr/bin/env python3

#read <dose> effect files and compute IC50

import sys, matplotlib.pylab as plt
import numpy as np
import math
from scipy.optimize import curve_fit

def sigmoid(x,top,bottom,slope,ic):
    return bottom + (top-bottom)/(1 + 10**((ic-x)*slope))

def sigmoid1(x,top,bottom,ic):
    return bottom + (top-bottom)/(1 + 10**((ic-x)))

data = np.loadtxt(sys.argv[1])
stage = 100
if len(sys.argv) > 2:
    stage = int(sys.argv[2])
    
x = data[:,0]
y = data[:,1]
ymax = np.max(y)
ymin = np.min(y)
mid = ymin + (ymax-ymin)/2.0
badline = np.polyfit(x,y,1)
badIC50 = (mid-badline[1])/badline[0]
print('Linear Fit IC50: %.2f'  % badIC50)


if stage >= 80:
    logx = np.log10(x)
    logline = np.polyfit(logx,y,1)
    logIC50 = (mid-logline[1])/logline[0]
    print('LogX-Linear Fit IC50: %.2f'  % 10**logIC50)

if stage >= 90:
    popt1,pcov1 = curve_fit(sigmoid1, logx, y)
    IC50_1 = popt1[2]
    print('Sigmoid Fit Fixed Slope IC50: %.2f' % 10**IC50_1)

if stage >= 100:
    popt,pcov = curve_fit(sigmoid, logx, y)
    IC50 = popt[3]
    slope = popt[2]
    print('Sigmoid Fit IC50: %.2f' % 10**IC50)
    print('Slope: %.3f' % slope)
