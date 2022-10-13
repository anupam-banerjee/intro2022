#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import sys

# 60: one file, plot only
# 70: one file, with labels and legend
# 80: two files
# 90: arbitrary number of files
# 95: arbitrary number of columns, show mean
# 100: arbitrary number of columns, show mean and stddev w/errorbars

# basic assignment:
# take three arguments, two input files and an output file name
# the input files have a title header and then x,y data
# output a linepoint graph

# bonus:
# take an arbitrary number of input files that have an arbitrary
# number of y value (e.g. x,y1,y2,y3)
# output a linepoint graph of the mean with the stddev as error bars
# unless there is only a single data point


stage = 100
try: #if last argument is number, assume it is specify stage
	stage = int(sys.argv[-1])
	sys.argv.pop()
except ValueError:
	pass
	
outname = sys.argv[-1]

end = -1
if stage < 80:
	end = 2
elif stage < 90:
	end = 3
	
for inname in sys.argv[1:end]:
	f = open(inname)
	header = f.readline().rstrip()
	data = np.array([x.rstrip().split(',') for x in f],float)
	xvals = data[:,0]
	
	if data.shape[1] > 2 and stage > 90: #do average and std w/errorbars
		yvals = np.mean(data[:,1:],axis=1)
		if stage > 95:
			dev = np.std(data[:,1:],axis=1)
			plt.errorbar(xvals,yvals,dev,fmt='-o',label=header)
		else:
			plt.plot(xvals,yvals,'-o',label=header)
	else:
		yvals = data[:,1]
		plt.plot(xvals,yvals,'-o',label=header)
		
if stage > 60:
	plt.xlabel("Time (min)")
	plt.ylabel("Cell Count")
	plt.legend(loc='best')
	
plt.savefig(outname)
