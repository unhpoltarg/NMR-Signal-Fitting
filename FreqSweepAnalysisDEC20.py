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

filename = 'Freq Sweeping\\simp120820.csv'
filename2 = 'Freq Sweeping\\simp120920.csv'
nmrwidth = 401
centFreq = 212.7
freqspan = 0.4
Toffset = 3690248400
#3690334786.06109
#3690334806.26108

file = open(filename, 'r')
reader = csv.reader(file)

baselines = []
TEdata = []
rawdata = []
times = []
datetime = []

rownum = 1
for row in reader:
    temp = []
    if rownum > 1:
        for i in range(nmrwidth):
            temp.append(float(row[6+i]))
        if rownum < 6:
            baselines.append(temp)
        if rownum > 5 and rownum < 770:
            TEdata.append(temp)
        if rownum > 769:
            datetime.append(row[0])
            times.append(float(row[1])-Toffset)
            rawdata.append(temp)
    rownum += 1

file.close()



file2 = open(filename2, 'r')
reader2 = csv.reader(file2)

rownum = 1
for row in reader2:
    temp = []
    if rownum > 1 and rownum < 1798:
        for i in range(nmrwidth):
            temp.append(float(row[10+i]))
        datetime.append(row[0])
        times.append(float(row[1])-Toffset)
        rawdata.append(temp)
    rownum += 1

file2.close()

freqs = []
baseavg = []
timesAvg = []

for i in range(nmrwidth):
    freqs.append(centFreq + freqspan*(i/(nmrwidth*1.0) - 0.5))

for i in range(nmrwidth):
    tempSum = 0.0
    for j in range(len(baselines)):
        tempSum += baselines[j][i]
    baseavg.append(tempSum/len(baselines))

newRaw = []
newTE = []
for j in range(len(rawdata)):
    temp = []
    for i in range(nmrwidth):
        temp.append(rawdata[j][i] - baseavg[i])
    newRaw.append(temp)

for j in range(len(TEdata)):
    temp = []
    for i in range(nmrwidth):
        temp.append(TEdata[j][i] - baseavg[i])
    newTE.append(temp)

wingFreqs = []
wingTE = []
wingData = []

for i in range(nmrwidth):
    if i < 70 or i > 330:
        wingFreqs.append(freqs[i])

for j in range(len(newTE)):
    temp = []
    for i in range(nmrwidth):
        if i < 70 or i > 330:
            temp.append(newTE[j][i])
    wingTE.append(temp)

for j in range(len(newRaw)):
    temp = []
    for i in range(nmrwidth):
        if i < 70 or i > 330:
            temp.append(newRaw[j][i])
    wingData.append(temp)

flatTE = []
flats = []

for j in range(len(newTE)):
    par, cov = curve_fit(deuts.cubic, wingFreqs, wingTE[j])
    temp = []
    for i in range(nmrwidth):
        temp.append(newTE[j][i] - deuts.cubic(freqs[i], par[0], par[1], par[2], par[3]))
    flatTE.append(temp)

for j in range(len(newRaw)):
    par, cov = curve_fit(deuts.cubic, wingFreqs, wingData[j])
    temp = []
    for i in range(nmrwidth):
        temp.append(newRaw[j][i] - deuts.cubic(freqs[i], par[0], par[1], par[2], par[3]))
    flats.append(temp)

TE_sums = []
TE_nums = []

flat_sums = []

for i in range(len(flatTE)):
    TE_sums.append(deuts.riemann_sum(freqs,flatTE[i]))
    TE_nums.append(i)

for i in range(len(flats)):
    flat_sums.append(deuts.riemann_sum(freqs,flats[i]))

TEcutoff = 200

totalTE = 0.0
TEcount = 0

for i in range(len(TE_sums)):
    if i > TEcutoff:
        totalTE += TE_sums[i]
        TEcount += 1

avgTE = totalTE/TEcount
print(avgTE)

ture = 1.1
magCur = 48.624
B = 0.1027405*magCur+0.00003497875
Pz_TE = deuts.TE_Pz_Proton(B,ture)
print(Pz_TE)

CC = Pz_TE/avgTE

Pzs = []
for i in range(len(flat_sums)):
    Pzs.append(flat_sums[i]*CC)

fig, az = plt.subplots()
az.scatter(times,Pzs,label='Enhanced Proton Pz, Dec 8/9 2020',color='black')
plt.xlabel("Seconds from Midnight 12/8/2020")
plt.ylabel("Pz")
plt.legend()
plt.show()

outstuff = []
for i in range(len(flats)):
    temp = []
    temp.append(datetime[i])
    temp.append(times[i])
    temp.append(Pzs[i])
    for k in range(nmrwidth):
        temp.append(flats[i][k])
    outstuff.append(temp)

with open('Freq Sweeping\\1220_Proton_flats.csv', mode='w', newline='', encoding='utf-8-sig') as Dec20FileP:
    writer = csv.writer(Dec20FileP, quoting=csv.QUOTE_ALL)
    writer.writerow(["DateTime","Sec from MN","Pz","Data"])
    writer.writerows(outstuff)
