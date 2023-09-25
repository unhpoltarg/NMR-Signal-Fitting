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

'''
Your data file should match ExampleCSV.csv:
Line 1: Central Frequency | Frequency Span | NMR Width (number of data points)
Line 2: Baseline data
Line 3: Signal data
YOU MUST HAVE THE deuts.py FILE IN THE SAME FOLDER
This is for fully real signals *only*.
'''

baseline = [] #A baseline signal taken when the magnet is not at larmor frequency
rawdata = [] #Your data

centfreq = 0.0
freqspan = 0.0
nmrwidth = 0

filename = 'ExampleCSV.csv' #Change to your csv's name - should be in same folder as this file
file = open(filename, 'r')
reader = csv.reader(file)

rownum = 1
for row in reader:
    if rownum == 1:
        centfreq = float(row[0])
        freqspan = float(row[1])
        nmrwidth = int(row[2])
    if rownum == 2:
        for i in range(nmrwidth):
            baseline.append(float(row[0+i]))
    if rownum == 3:
        for i in range(nmrwidth):
            rawdata.append(float(row[0+i]))
    rownum += 1

file.close()

freqs = [] #This propogates the frequencies associated with your data points

for i in range(nmrwidth):
    freqs.append(centfreq + freqspan*(i/(nmrwidth*1.0) - 0.5))

newRaw = [] #Subtracting the baseline from the signal

for i in range(nmrwidth):
    newRaw.append(rawdata[i] - baseline[i])

wingFreqs = [] #Fitting the edges of the signal to flatten the sides
wingData = []

for i in range(nmrwidth):
    if i < nmrwidth*0.2 or i > nmrwidth*0.8:
        wingFreqs.append(freqs[i])
            
for j in range(len(newRaw)):
    temp = []
    for i in range(nmrwidth):
        if i < nmrwidth*0.2 or i > nmrwidth*0.8:
            temp.append(newRaw[i])
    wingData.append(temp)

flats = []

par, cov = curve_fit(deuts.cubic, wingFreqs, wingData[j])
for i in range(nmrwidth):
    flats.append(newRaw[i] - deuts.cubic(freqs[i], par[0], par[1], par[2], par[3]))

#This will produce a plot of your flattened signal. If the sides are not flat, or the signal
# looks wrong, go back and check your csv to make sure it adheres to requirements.
fig, ax = plt.subplots()
ax.scatter(freqs,flats,label='Flattened Signal',color='black')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

#This fits the flattened signal. The parameters will print out in IDLE.
#If any parameters are "railing" (going to X.99999 or X.00001), adjust the minima or maxima.
pFit, cFit = curve_fit(deuts.deut_double, freqs, flats, maxfev=100000,
                           bounds=([centfreq-0.15,0.0214,0.025,0.01,0,0.0001,1.1,0,0],
                               [centfreq+0.15,0.025,0.035,0.06,0.1,0.4,2.5,0.4,0.005]))
temp = []
temp.append(pFit[0]) # Omega_D - the central frequency of the signal
# Omega_Q - the quadrupole frequency: 3*omegaQ = distance from omegaD to the peak
temp.append(pFit[1]) # Omega_Q1 - carbon quadrupole frequency
temp.append(pFit[2]) # Omega_Q2 - oxygen quadrupole frequency
# A - the Lorentzian width of the signal
temp.append(pFit[3]) # A1 - carbon Lorentzian width (oxygen's is pegged to it: omegaQ1*A1 = omegaQ2*A2)
# eta - asymmetry parameter
temp.append(pFit[4]) # eta1 - carbon
temp.append(pFit[5]) # eta2 - oxygen
temp.append(pFit[6]) # r - "ratio" number that describes Pz and Pzz
temp.append(pFit[7]) # K - the percentage of deuteron bonds that are with oxygen
temp.append(pFit[8]) # scaling factor - turns the proportionality of chi'' into a function
print(temp)

#designate the fit
fit = []
for i in range(nmrwidth):
    fit.append(deuts.deut_double(freqs[i],pFit[0],pFit[1],pFit[2],pFit[3],pFit[4],pFit[5],pFit[6],pFit[7],pFit[8]))

#Print in IDLE Pz and Pzz in decimal form.
r = pFit[6]
Pz = deuts.ratio_method(r)
Pzz = deuts.ratio_tensor(r)
print("Pz = " + str(Pz))
print("Pzz = " + str(Pzz))

#This shows a plot of your fit compared to the data.
fig, ax = plt.subplots()
ax.scatter(freqs,flats,label='Data',color='black')
ax.plot(freqs,fit,label='Fit',color='blue')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show() 
