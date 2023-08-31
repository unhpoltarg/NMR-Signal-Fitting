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

filename = 'Freq Sweeping\\simp050622REAL.csv'
#filename2 = 'Freq Sweeping\\simp91520.csv'
nmrwidth = 451
#centFreq = 33.05
freqspan = 0.8
Toffset = 0

file = open(filename, 'r')
reader = csv.reader(file)

baselines = []
rawdata = []
times = []
datetime = []
CF = 32.792

rownum = 1
for row in reader:
    temp = []
    if rownum > 1:
        for i in range(nmrwidth):
            temp.append(float(row[8+i]))
        if rownum == 196 or rownum == 197 or rownum == 198:
            baselines.append(temp)
        if rownum < 196 or rownum > 198 and rownum < 625:
            datetime.append(row[0])
            times.append(float(row[1])-Toffset)
            rawdata.append(temp)
    rownum += 1

file.close()

freqs = []
baseavg = []
dataAvg = []
timesAvg = []


for i in range(nmrwidth):
    freqs.append(CF + freqspan*(i/(nmrwidth*1.0) - 0.5))
    tempSum = 0.0
    for k in range(len(baselines)):
        tempSum += baselines[k][i]
    baseavg.append(tempSum/len(baselines))
'''
fig, ax = plt.subplots()
ax.scatter(freqs,baselines[0],label='0',color='black')
ax.scatter(freqs,baselines[1],label='1',color='red')
ax.scatter(freqs,baselines[2],label='2',color='orange')
ax.scatter(freqs,baselines[3],label='3',color='yellow')
ax.scatter(freqs,baselines[4],label='4',color='green')
plt.legend()
plt.show()
'''
newRaw = []

for j in range(len(rawdata)):
    temp = []
    for i in range(nmrwidth):
        temp.append(rawdata[j][i] - baseavg[i])
    newRaw.append(temp)

wingFreqs = []
wingData = []

for i in range(nmrwidth):
    if i < 70 or i > 381:
        wingFreqs.append(freqs[i])
            
for j in range(len(newRaw)):
    temp = []
    for i in range(nmrwidth):
        if i < 70 or i > 381:
            temp.append(newRaw[j][i])
    wingData.append(temp)

flats = []

for j in range(len(newRaw)):
    par, cov = curve_fit(deuts.cubic, wingFreqs, wingData[j])
    temp = []
    for i in range(nmrwidth):
        temp.append(newRaw[j][i] - deuts.cubic(freqs[i], par[0], par[1], par[2], par[3]))
    flats.append(temp)
'''
print(times[321])
fig, ax = plt.subplots()
ax.scatter(freqs,flats[100],label='100',color='black')
ax.scatter(freqs,flats[200],label='200',color='red')
ax.scatter(freqs,flats[400],label='400',color='orange')
ax.scatter(freqs,flats[600],label='600',color='yellow')
ax.hlines(y=0,xmin=freqs[0],xmax=freqs[len(freqs)-1],color='blue')
plt.legend()
plt.show()
'''
'''
pFit, cFit = curve_fit(deuts.deut_double, freqs, flats[600], maxfev=100000,
                           bounds=([32.75,0.0214,0.025,0.01,0,0.0001,1.1,0,0],
                               [33.0,0.025,0.035,0.06,0.1,0.4,2.0,0.4,0.005]))
temp = []
temp.append(pFit[0])
temp.append(pFit[1])
temp.append(pFit[2])
temp.append(pFit[3])
temp.append(pFit[4])
temp.append(pFit[5])
temp.append(pFit[6])
temp.append(pFit[7])
temp.append(pFit[8])

temp2 = []
print(temp)

for i in range(nmrwidth):
    temp2.append(deuts.deut_double(freqs[i],pFit[0],pFit[1],pFit[2],pFit[3],pFit[4],pFit[5],pFit[6],pFit[7],pFit[8]))

fig, ax = plt.subplots()
ax.scatter(freqs,flats[600],label='Data',color='black')
ax.plot(freqs,temp2,label='Fit',color='blue')
plt.legend()
plt.show() 
'''

fits = []
rS = []

for j in range(len(flats)):
    print(str(j+1) + " out of " + str(len(flats)+1))
    temp = []
    pS, cS = curve_fit(deuts.DD_hardcode, freqs, flats[j], maxfev=100000,
                    bounds=([0.0],[2.0]))
    for i in range(nmrwidth):
        temp.append(deuts.DD_hardcode(freqs[i],pS[0]))
    fits.append(temp)
    print("r = " + str(pS[0]))
    rS.append(pS[0])
  
Pz = []
Pzz = []

for i in range(len(rS)):
    Pz.append(deuts.ratio_method(rS[i]))
    Pzz.append(deuts.ratio_tensor(rS[i]))

fig, ax = plt.subplots()
ax.scatter(times,Pz,label='Pz',color='black')
plt.xlabel("Seconds from Midnight 5/6/2022")
plt.ylabel("Pz")
plt.legend()
plt.show()

outstuff = []
for i in range(len(fits)):
    temp = []
    temp.append(datetime[i])
    temp.append(times[i])
    temp.append(rS[i])
    temp.append(Pz[i])
    temp.append(Pzz[i])
    for k in range(nmrwidth):
        temp.append(fits[i][k])
    outstuff.append(temp)

with open('Freq Sweeping\\050622fits.csv', mode='w', newline='', encoding='utf-8-sig') as May22File:
    writer = csv.writer(May22File, quoting=csv.QUOTE_ALL)
    writer.writerow(["DateTime","Sec from MN","r","Pz","Pzz","Data"])
    writer.writerows(outstuff)

outstuff2 = []
for i in range(len(flats)):
    temp = []
    temp.append(datetime[i])
    temp.append(times[i])
    temp.append(rS[i])
    temp.append(Pz[i])
    temp.append(Pzz[i])
    for k in range(nmrwidth):
        temp.append(flats[i][k])
    outstuff2.append(temp)

with open('Freq Sweeping\\050622flats.csv', mode='w', newline='', encoding='utf-8-sig') as May22File2:
    writer2 = csv.writer(May22File2, quoting=csv.QUOTE_ALL)
    writer2.writerow(["DateTime","Sec from MN","r","Pz","Pzz","Data"])
    writer2.writerows(outstuff2)
