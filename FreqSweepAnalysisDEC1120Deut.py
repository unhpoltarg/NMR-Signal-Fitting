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

filename = 'Freq Sweeping\\simp121020.csv'
filename2 = 'Freq Sweeping\\simp121120.csv'
nmrwidth = 401
centFreq = 33.05
freqspan = 0.8
Toffset = 3690507600


file = open(filename, 'r')
reader = csv.reader(file)

baselines = []
rawdata = []
times = []
datetime = []

rownum = 1
for row in reader:
    temp = []
    if rownum > 1:
        for i in range(nmrwidth):
            temp.append(float(row[10+i]))
        if rownum < 11:
            baselines.append(temp)
    rownum += 1

file.close()

file2 = open(filename2, 'r')
reader2 = csv.reader(file2)

rownum = 1
for row in reader2:
    temp = []
    if rownum > 1:
        datetime.append(row[0])
        times.append(float(row[1])-Toffset)
        for i in range(nmrwidth):
            temp.append(float(row[10+i]))
        rawdata.append(temp)
    rownum += 1

file2.close()

freqs = []
baseavg = []
dataAvg = []
timesAvg = []

for i in range(nmrwidth):
    freqs.append(centFreq + freqspan*(i/(nmrwidth*1.0) - 0.5))
    tempSum = 0.0
    for j in range(len(baselines)):
        tempSum += baselines[j][i]
    baseavg.append(tempSum/len(baselines))


newRaw = []
for j in range(len(rawdata)):
    temp = []
    for i in range(nmrwidth):
        temp.append(rawdata[j][i] - baseavg[i])
    newRaw.append(temp)

wingFreqs = []
wingData = []

for i in range(nmrwidth):
    if i < 70 or i > 330:
        wingFreqs.append(freqs[i])

for j in range(len(newRaw)):
    temp = []
    for i in range(nmrwidth):
        if i < 70 or i > 330:
            temp.append(newRaw[j][i])
    wingData.append(temp)

flats = []

for j in range(len(newRaw)):
    par, cov = curve_fit(deuts.cubic, wingFreqs, wingData[j])
    temp = []
    for i in range(nmrwidth):
        temp.append(newRaw[j][i] - deuts.cubic(freqs[i], par[0], par[1], par[2], par[3]))
    flats.append(temp)

'''fig, ax = plt.subplots()
ax.scatter(freqs,flats[1700],label='0',color='black')
ax.scatter(freqs,flats[100],label='1',color='red')
ax.scatter(freqs,flats[200],label='2',color='orange')
ax.scatter(freqs,flats[300],label='3',color='yellow')
ax.scatter(freqs,flats[400],label='4',color='green')
ax.scatter(freqs,flats[500],label='5',color='cyan')
ax.scatter(freqs,flats[600],label='6',color='blue')
ax.scatter(freqs,flats[700],label='7',color='purple')
ax.scatter(freqs,flats[800],label='8',color='pink')
plt.legend()
plt.show()'''

'''
pFit, cFit = curve_fit(deuts.deut_double, freqs, flats[1700], maxfev=100000,
                           bounds=([32.7,0.0214,0.026,0.01,0,0,1.1,0,0],
                               [33.3,0.022,0.029,0.04,0.1,0.2,3.0,0.2,0.05]))
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
ax.scatter(freqs,flats[1700],label='Data',color='black')
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
plt.xlabel("Seconds from Midnight 12/11/2020")
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

with open('Freq Sweeping\\121120fits.csv', mode='w', newline='', encoding='utf-8-sig') as Dec20File:
    writer = csv.writer(Dec20File, quoting=csv.QUOTE_ALL)
    writer.writerow(["DateTime","Sec from MN","r","Pz","Data"])
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

with open('Freq Sweeping\\121120flats.csv', mode='w', newline='', encoding='utf-8-sig') as Dec20File2:
    writer = csv.writer(Dec20File2, quoting=csv.QUOTE_ALL)
    writer.writerow(["DateTime","Sec from MN","r","Pz","Pzz","Data"])
    writer.writerows(outstuff2)
