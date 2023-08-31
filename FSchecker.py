import math
import plotly
import plotly.graph_objs as go
import numpy as np
import matplotlib.pyplot as plt
import csv
from datetime import datetime
from scipy.optimize import curve_fit
from datetime import datetime
import deuts

filename = 'Freq Sweeping\\920flats.csv'
filename2 = 'Freq Sweeping\\920fits.csv'
nmrwidth = 401
centFreq = 32.925
freqspan = 0.8

flats = []
fits = []
datetime = []
time = []

file = open(filename, 'r')
reader = csv.reader(file)

rownum = 1
for row in reader:
    temp = []
    if rownum > 1:
        for i in range(nmrwidth):
            temp.append(float(row[4+i]))
        datetime.append(row[0])
        time.append(row[1])
        flats.append(temp)
    rownum += 1

file.close()

file2 = open(filename2, 'r')
reader2 = csv.reader(file2)

rownum = 1
for row in reader2:
    temp = []
    if rownum > 1:
        for i in range(nmrwidth):
            temp.append(float(row[4+i]))
        fits.append(temp)
    rownum += 1

file.close()

freqs = []
for i in range(nmrwidth):
    freqs.append(centFreq + freqspan*(i/(nmrwidth*1.0) - 0.5))

fig, ax = plt.subplots()
ax.scatter(freqs,flats[3703],label='Data',color='black')
ax.plot(freqs,fits[3703],label='Fit',color='red',linestyle='solid')
plt.legend()
plt.title("Signal from " + str(datetime[3703]))
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.show()

fig, ax = plt.subplots()
ax.scatter(freqs,flats[3704],label='Data',color='black')
ax.plot(freqs,fits[3704],label='Fit',color='red',linestyle='solid')
plt.legend()
plt.title("Signal from " + str(datetime[3704]))
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.show() 
