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

filename = 'Freq Sweeping\\simp100121.csv'
#filename2 = 'Freq Sweeping\\simp91520.csv'
nmrwidth = 401
#centFreq = 33.05
freqspan = 0.8
Toffset = 0

file = open(filename, 'r')
reader = csv.reader(file)

baselines = []
rawdata = []
times = []
datetime = []
CFs = [33.05,33.626,33.56,33]

rownum = 1
for row in reader:
    temp = []
    if rownum > 1:
        for i in range(nmrwidth):
            temp.append(float(row[6+i]))
        if rownum < 4:
            baselines.append(temp)
        if rownum > 3 and rownum < 648:
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
    freqs.append(CFs[0] + freqspan*(i/(nmrwidth*1.0) - 0.5))
    tempSum = 0.0
    for k in range(len(baselines)):
        tempSum += baselines[k][i]
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
    if i < 70 or i > 345:
        wingFreqs.append(freqs[i])
            
for j in range(len(newRaw)):
    temp = []
    for i in range(nmrwidth):
        if i < 70 or i > 345:
            temp.append(newRaw[j][i])
    wingData.append(temp)

flats = []

for j in range(len(newRaw)):
    par, cov = curve_fit(deuts.cubic, wingFreqs, wingData[j])
    temp = []
    for i in range(nmrwidth):
        temp.append(newRaw[j][i] - deuts.cubic(freqs[i], par[0], par[1], par[2], par[3]))
    flats.append(temp)

print(times[321])
fig, ax = plt.subplots()
ax.scatter(freqs,rawdata[321],label='321',color='black')
ax.scatter(freqs,flats[322],label='322',color='red')
ax.scatter(freqs,flats[323],label='323',color='orange')
ax.scatter(freqs,flats[324],label='324',color='yellow')
ax.scatter(wingFreqs,wingData[320],label='320 wings',color='green')
plt.legend()
plt.show()

'''
pFit, cFit = curve_fit(deuts.deut_double, freqs[0], flats[0][300], maxfev=100000,
                           bounds=([33.01,0.0214,0.025,0.01,0,0.0001,1.1,0,0],
                               [33.3,0.025,0.035,0.06,0.1,0.4,3.0,0.4,0.05]))
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
    temp2.append(deuts.deut_double(freqs[0][i],pFit[0],pFit[1],pFit[2],pFit[3],pFit[4],pFit[5],pFit[6],pFit[7],pFit[8]))

fig, ax = plt.subplots()
ax.scatter(freqs[0],flats[0][300],label='Data',color='black')
ax.plot(freqs[0],temp2,label='Fit',color='blue')
plt.legend()
plt.show() 
'''
'''
fits = [[] for i in range(4)]
rS = [[] for i in range(4)]

for k in range(4):
    for j in range(len(flats[k])):
        print(str(j+1) + " out of " + str(len(flats[k])+1))
        temp = []
        pS, cS = curve_fit(deuts.DD_HC_OpenD, freqs[k], flats[k][j], maxfev=100000,
                       bounds=([32.8,0.0],[33.8,2.0]))
        for i in range(nmrwidth):
            temp.append(deuts.DD_HC_OpenD(freqs[k][i],pS[0],pS[1]))
        fits[k].append(temp)
        print("r = " + str(pS[1]))
        rS[k].append(pS[1])
  
Pz = [[] for i in range(4)]
Pzz = [[] for i in range(4)]

for k in range(4):
    for i in range(len(rS[k])):
        Pz[k].append(deuts.ratio_method(rS[k][i]))
        Pzz[k].append(deuts.ratio_tensor(rS[k][i]))

fig, ax = plt.subplots()
ax.scatter(times[0],Pz[0],label='Pz0',color='black')
ax.scatter(times[1],Pz[1],label='Pz1',color='blue')
ax.scatter(times[2],Pz[2],label='Pz2',color='green')
ax.scatter(times[3],Pz[3],label='Pz3',color='red')
plt.xlabel("Seconds from Midnight 12/8/2020")
plt.ylabel("Pz")
plt.legend()
plt.show()

outstuff = []
for j in range(4):
    for i in range(len(fits[j])):
        temp = []
        temp.append(datetime[j][i])
        temp.append(times[j][i])
        temp.append(CFs[j])
        temp.append(rS[j][i])
        temp.append(Pz[j][i])
        temp.append(Pzz[j][i])
        for k in range(nmrwidth):
            temp.append(fits[j][i][k])
        outstuff.append(temp)

with open('Freq Sweeping\\100121fits.csv', mode='w', newline='', encoding='utf-8-sig') as Oct21File:
    writer = csv.writer(Oct21File, quoting=csv.QUOTE_ALL)
    writer.writerow(["DateTime","Sec from MN","Central Freq","r","Pz","Pzz","Data"])
    writer.writerows(outstuff)

outstuff2 = []
for j in range(4):
    for i in range(len(flats[j])):
        temp = []
        temp.append(datetime[j][i])
        temp.append(times[j][i])
        temp.append(CFs[j])
        temp.append(rS[j][i])
        temp.append(Pz[j][i])
        temp.append(Pzz[j][i])
        for k in range(nmrwidth):
            temp.append(flats[j][i][k])
        outstuff2.append(temp)

with open('Freq Sweeping\\100121flats.csv', mode='w', newline='', encoding='utf-8-sig') as Oct21File2:
    writer2 = csv.writer(Oct21File2, quoting=csv.QUOTE_ALL)
    writer2.writerow(["DateTime","Sec from MN","Central Freq","r","Pz","Pzz","Data"])
    writer2.writerows(outstuff2)
'''
