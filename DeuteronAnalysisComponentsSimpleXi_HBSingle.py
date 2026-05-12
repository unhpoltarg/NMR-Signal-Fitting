import math
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
***
'''
#YOU MUST HAVE THE deuts.py FILE IN THE SAME FOLDER

#expected values for propanediol
'''omegQ1 = 2.17139774e-02
omegQ2 = 2.88692494e-02
A1 = 3.06372826e-02
eta1 = 4.59723046e-02
eta2 = 6.12644049e-02
K = 1.15523343e-01'''

#expected values for TEMPO butanol
'''omegQ1 = 2.17004538e-02
omegQ2 = 2.89859414e-02
A1 = 2.75141728e-02
eta1 = 6.72464295e-02
eta2 = 6.57934478e-02
K = 4.28479420e-02'''

#expected values for irradiated butanol
omegQ1 = 0.027335
A1 = 0.038029
eta1 = 0.116228
scale = 0.001189
#K = 0.136 according to Dulya, seems wrong
#scale = 9.66383779e-04

'''def DRDC_HC(freq,omegD,r,scale,phase,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase,x3,x2,x1,x0)
    return out

def DRDC_Xi_HC(freq,omegD,r,scale,phase,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase,x3,x2,x1,x0)
    return out'''

def DRDC_XiO_HC(freq,omegD,r,phaseX,x3,x2,x1,x0):
    out = deuts.deut_real_cube_Xi(freq,omegD,omegQ1,A1,eta1,r,scale,phaseX,x3,x2,x1,x0)
    return out

def DRDC_XiO_scHC(freq,omegD,r,scaleX,phaseX,x3,x2,x1,x0):
    out = deuts.deut_real_cube_Xi(freq,omegD,omegQ1,A1,eta1,r,scaleX,phaseX,x3,x2,x1,x0)
    return out

def DRD_XiO_scHC(freq,omegD,r,scaleX,phaseX):
    out = deuts.deut_real_Xi(freq,omegD,omegQ1,A1,eta1,r,scaleX,phaseX)
    return out

print(str(deuts.ratio_method(2.5)))
bases = []
raws = []
HBs = []
rawdata = [] #Your equilibrium data
HBdata = [] #Your holeburn data
baseline = [] #your baseline
holeline = []

centfreq = 0.0
freqspan = 0.0
nmrwidth = 0

filename = 'ExampleCSV_HBavg.csv' #Change to your csv's name - should be in same folder as this file
file = open(filename, 'r')
reader = csv.reader(file)

rownum = 1
for row in reader:
    if rownum == 1:
        centfreq = float(row[0])
        freqspan = float(row[1])
        nmrwidth = int(row[2])
    if rownum >= 2 and rownum <= 16:
        temp = []
        for i in range(nmrwidth):
            temp.append(float(row[0+i]))
        bases.append(temp)
    if rownum >= 18 and rownum <= 21:
        temp = []
        for i in range(nmrwidth):
            temp.append(float(row[0+i]))
        raws.append(temp)
    if rownum >= 23 and rownum <= 25:
        temp = []
        for i in range(nmrwidth):
            temp.append(float(row[0+i]))
        HBs.append(temp)
    '''if rownum == 5:
        for i in range(nmrwidth):
            holeline.append(float(row[0+i]))'''
    rownum += 1

file.close()

print("Cent freq = " + str(centfreq))

freqs = [] #This propogates the frequencies associated with your data points

for i in range(nmrwidth):
    freqs.append(centfreq + freqspan*(i/(nmrwidth*1.0) - 0.5))
    temp = 0.0
    for j in range(len(bases)):
        temp += bases[j][i]
    baseline.append(temp/len(bases))
    temp = 0.0
    for j in range(len(raws)):
        temp += raws[j][i]
    rawdata.append(temp/len(raws))
    temp = 0.0
    for j in range(len(HBs)):
        temp += HBs[j][i]
    HBdata.append(temp/len(HBs))

newRaw = [] #Subtracting the baseline from the signals
newHB = []

for i in range(nmrwidth):
    newRaw.append(rawdata[i] - baseline[i])
    newHB.append(HBdata[i] - baseline[i])# - holeline[i])#

freqs = freqs[80:nmrwidth - 80]
newRaw = newRaw[80:nmrwidth - 80]
newHB = newHB[80:nmrwidth - 80]
nmrwidth = len(freqs)

#This will produce a plot of your flattened signal. If the signal
# looks wrong, go back and check your csv to make sure it adheres to requirements.
fig, ax = plt.subplots()
ax.scatter(freqs,newRaw,label='Eq Signal',color='black')
ax.scatter(freqs,newHB,label='HB Signal',color='red')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

pFit2, cFit2 = curve_fit(DRDC_XiO_HC, freqs, newRaw, maxfev=100000,
                           bounds=([centfreq-0.5,1,-2,-1000,-1000,-10000,-100000],
                               [centfreq+0.5,2,2,1000,1000,10000,100000]))

print(pFit2)
omegD = pFit2[0]
rX = pFit2[1]
#scaleX = pFit2[2]
phase = pFit2[2]

x3x = pFit2[3]
x2x = pFit2[4]
x1x = pFit2[5]
x0x = pFit2[6]

cube_subtX = []
fitX = []
for i in range(nmrwidth):
    cube_subtX.append(newRaw[i] - deuts.cubic(freqs[i],x3x,x2x,x1x,x0x))
    fitX.append(deuts.deut_real_Xi(freqs[i],omegD,omegQ1,A1,eta1,rX,scale,phase))

fig, ax = plt.subplots()
ax.scatter(freqs,cube_subtX,label='Eq Signal',color='black')
ax.plot(freqs,fitX,label='Fit',color='red')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

#wing fit
wingF = []
wingD = []
for i in range(nmrwidth):
    if i < 30 or i > nmrwidth - 30:
        wingF.append(freqs[i])
        wingD.append(newHB[i])
pWing, cWing = curve_fit(deuts.cubic, wingF, wingD, maxfev=100000,
                         bounds=([-1000,-1000,-10000,-100000],
                               [1000,1000,10000,100000]))

x3HB = pWing[0]
x2HB = pWing[1]
x1HB = pWing[2]
x0HB = pWing[3]

'''pFit3, cFit3 = curve_fit(DRD_XiO_scHC, freqs, newHB, maxfev=100000,
                           bounds=([centfreq-0.5,1,0,-2],
                               [centfreq+0.5,20,1,2]))

print(pFit3)
omegDHB = pFit3[0]
rHB = pFit3[1]
scaleHB = pFit3[2]
phaseHB = pFit3[3]'''

'''#HB fit
pFit3, cFit3 = curve_fit(DRDC_XiO_scHC, freqs, newHB, maxfev=100000,
                           bounds=([centfreq-0.5,1,0,-2,-1000,-1000,-10000,-100000],
                               [centfreq+0.5,20,1,2,1000,1000,10000,100000]))

print(pFit3)
omegDHB = pFit3[0]
rHB = pFit3[1]
scaleHB = pFit3[2]
phaseHB = pFit3[3]

x3HB = pFit3[4]
x2HB = pFit3[5]
x1HB = pFit3[6]
x0HB = pFit3[7]'''

cube_subtHB = []
#cube_subtHBa = []
prefit = []
fitHB = []
diff = []

for i in range(nmrwidth):
    cube_subtHB.append(newHB[i] - deuts.cubic(freqs[i],x3HB,x2HB,x1HB,x0HB))
    #fitHB.append(deuts.deut_real_double_Xi(freqs[i],omegDHB,omegQ1,omegQ2,A1,eta1,eta2,rHB,K,scaleHB,phaseHB))
    #prefit.append(DRDC_XiO_scHC(freqs[i],omegDHB,rHB,scaleHB,phaseHB,x3HB,x2HB,x1HB,x0HB))
    
#scaleNu = cube_subtX[226]/cube_subtHBa[226]
for i in range(nmrwidth):
    #cube_subtHB.append(cube_subtHBa[i]*scaleNu)
    diff.append(cube_subtHB[i] - cube_subtX[i])

'''fig, ax = plt.subplots()
ax.scatter(freqs,newHB,label='HB Signal',color='black')
ax.plot(freqs,prefit,label='Pre-Fit',color='red')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()'''

fig, ax = plt.subplots()
ax.scatter(freqs,cube_subtHB,label='HB Signal',color='black')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

#look at cubic subtracted differences
fig, ax = plt.subplots()
#ax.scatter(freqs,cube_subtX,label='Eq Signal',color='black')
#ax.scatter(freqs,cube_subtHB,label='HB Signal',color='red')
ax.scatter(freqs,diff,label='Difference',color='blue')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

rZ = []
rZ.append(rX)
rBump = []
rBump.append([1.0])

def HB_HC(freq,arr):
    outt = deuts.deut_real_Xi(freq,omegD,omegQ1,A1,eta1,arr,scale,phase)
    return outt

fitHB = []
fitHB_IM = []
fitComps = [[] for i in range(4)]
fitCompsIm = [[] for i in range(4)]
eqComps = [[] for i in range(4)]
eqCompsIm = [[] for i in range(4)]

def Bump(freq,req,arr):
            outt = (deuts.deut_real_Xi(freq,omegD,omegQ1,A1,eta1,req,scale,phase) +
                    0.5*deuts.deut_real_Xi_rev(freq,omegD,omegQ1,A1,eta1,req,scale,phase) -
                    0.5*deuts.deut_real_Xi_rev(freq,omegD,omegQ1,A1,eta1,arr,scale,phase))                                                                                                        
            return outt

def NormComp(freq,arr,phi,epsi):
    omegQ = omegQ1
    A = A1
    eta = eta1
    out = deuts.Xi(arr)*scale*(deuts.deut_phiavg(freq,epsi,omegD,omegQ,A,eta,arr)*np.cos(phi) + 
            deuts.deut_disp_phiavg(freq,epsi,omegD,omegQ,A,eta,arr)*np.sin(phi))
    return out

def NormCompIm(freq,arr,epsi):
    omegQ = omegQ1
    A = A1
    eta = eta1
    out = deuts.Xi(arr)*scale*(deuts.deut_phiavg(freq,epsi,omegD,omegQ,A,eta,arr)*np.cos(math.pi/2.0) + 
            deuts.deut_disp_phiavg(freq,epsi,omegD,omegQ,A,eta,arr)*np.sin(math.pi/2.0))
    return out

def BumpComp(freq,req,rstc,phi,epsi): #CO = C or O
    omegQ = omegQ1
    A = A1
    eta = eta1
    equi = NormComp(freq,req,phi,epsi)
    fake = 0.5*deuts.Xi(req)*scale*(deuts.deut_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,req)*np.cos(phi) + 
                                deuts.deut_disp_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,req)*np.sin(phi))
    holeSub = 0.5*deuts.Xi(rstc)*scale*(deuts.deut_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,rstc)*np.cos(phi) + 
                                deuts.deut_disp_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,rstc)*np.sin(phi))
    out = equi + fake - holeSub
    return out

def BumpCompIm(freq,req,rstc,epsi): #CO = C or O
    omegQ = omegQ1
    A = A1
    eta = eta1
    equi = NormCompIm(freq,req,epsi)
    fake = 0.5*deuts.Xi(req)*scale*(deuts.deut_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,req)*np.cos(math.pi/2.0) + 
                                deuts.deut_disp_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,req)*np.sin(math.pi/2.0))
    holeSub = 0.5*deuts.Xi(rstc)*scale*(deuts.deut_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,rstc)*np.cos(math.pi/2.0) + 
                                deuts.deut_disp_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,rstc)*np.sin(math.pi/2.0))
    out = equi + fake - holeSub
    return out

def BumpHC(freq,arr):
    out = Bump(freq,rX,arr)
    return out
for i in range(nmrwidth-2):
    print(str(i+1) + "/" + str(nmrwidth-2))
    tempD = []
    tempD.append(cube_subtHB[0])
    tempD.append(cube_subtHB[i+1])
    tempD.append(cube_subtHB[nmrwidth-1])
    tempF = []
    tempF.append(freqs[0])
    tempF.append(freqs[i+1])
    tempF.append(freqs[nmrwidth-1])
    pHB, cHB = curve_fit(HB_HC, tempF, tempD, maxfev=10000,
                           bounds=([1],
                               [3.5]))
    print(pHB[0])
    rZ.append(pHB[0])
        
    if rZ[i] > rX:
        print("Bump")
        pB, cB = curve_fit(BumpHC, tempF, tempD, maxfev=10000,
                           bounds=([1],
                               [3]))
        rBump.append(pB[0])
    else:
        rBump.append(pHB[0])
rZ.append(rX)

fig, ax = plt.subplots()
ax.scatter(freqs,rZ,label='r',color='green')
plt.xlabel("Frequency (MHz)")
plt.ylabel("P")
plt.legend()
plt.show()

rBump.append(1.0)
fitHB = []
for i in range(len(freqs)):
        fitHB.append(deuts.deut_real_Xi(freqs[i],omegD,omegQ1,A1,eta1,rZ[i],scale,phase))
        #fitHB_IM.append(deuts.deut_real_double_Xi(freqs[i],omegD,omegQ1,omegQ2,A1,eta1,eta2,rZ[i],K,scaleX,math.pi/2.0))

for i in range(nmrwidth):
    eqComps[0].append(NormComp(freqs[i],rX,phase,-1))
    eqComps[1].append(NormComp(freqs[i],rX,phase,1))
    if rZ[i] > rX:
        fitComps[0].append(BumpComp(freqs[i],rX,rBump[i],phase,-1))
        fitComps[1].append(BumpComp(freqs[i],rX,rBump[i],phase,1))           
    else:
        fitComps[0].append(NormComp(freqs[i],rZ[i],phase,-1))
        fitComps[1].append(NormComp(freqs[i],rZ[i],phase,1))

#simEq = []
fitHB2 = []
for i in range(len(rZ)):
    fitHB2.append(deuts.deut_real_Xi(freqs[i],omegD,omegQ1,A1,eta1,rZ[i],scale,phase))
    #simEq.append(deuts.deut_real_double_Xi(freqs[i],omegDX,omegQ1,omegQ2,A1,eta1,eta2,rZ[226],K,scaleX,phaseHB))

#Print in IDLE Pz and Pzz in decimal form.
STD = 0.0009863614666752998

PzX = deuts.ratio_method(rX)
PzzX = deuts.ratio_tensor(rX)
print("P = " + str(100*PzX)+"%")
print("Q = " + str(100*PzzX)+"%")
CCeq = PzX/deuts.riemann_sum(freqs,cube_subtX)
eqSTD = 0.0
for i in range(nmrwidth):
    eqSTD += abs(fitX[i] - cube_subtX[i])
eqSTD = (eqSTD/nmrwidth)*CCeq*100
print("Eq STD = " + str(eqSTD))

print("Error = " + str(eqSTD)+"%")
print("Phase angle = " + str((180*phase)/np.pi) + " degrees")
print("Reduced chi-squared = " + str(deuts.chisq_red(cube_subtX,fitX,STD)))
print("CC = " + str(CCeq))

IminC = deuts.riemann_sum(freqs,fitComps[0])
IplusC = deuts.riemann_sum(freqs,fitComps[1])
HB_P = CCeq*(IplusC + IminC)
HB_Q = CCeq*(IplusC - IminC)
    
'''holes = 0.0
bumps = 0.0
for i in range(nmrwidth):
        if rZ[i] < rEQX:
            holes += eqComps[j][0][i] - fitComps[j][0][i]
            holes += eqComps[j][1][i] - fitComps[j][1][i]
            holes += eqComps[j][2][i] - fitComps[j][2][i]
            holes += eqComps[j][3][i] - fitComps[j][3][i]
        if rZ[j][i] > rEQ[j]:
            bumps += fitComps[j][0][i] - eqComps[j][0][i]
            bumps += fitComps[j][1][i] - eqComps[j][1][i]
            bumps += fitComps[j][2][i] - eqComps[j][2][i]
            bumps += fitComps[j][3][i] - eqComps[j][3][i]'''

print("Hole burnt P = " +str(100*HB_P)+"%")
print("Hole burnt Q = " +str(100*HB_Q)+"%")

'''    print("Holes total = " +str(holes))
    print("Bumps total = " +str(bumps))'''

diff2 = []
for i in range(nmrwidth):
    diff2.append(fitHB2[i] - fitX[i])

'''added = []
addedQ = []
for i in range(nmrwidth):
    temp = 0.0
    for j in range(len(fitComps)):
        temp += fitComps[i]
    added.append(temp)
    temp = 0.0
    for j in range(len(fitCompsQ)):
        temp += fitCompsQ[i]
    addedQ.append(temp)'''



'''with open('Feb13_2025_HB.csv', mode='w', newline='', encoding='utf-8-sig') as May22File:
    writer = csv.writer(May22File, quoting=csv.QUOTE_ALL)
    writer.writerow(["HB Data"])
    writer.writerow(cube_subtHB)
    writer.writerow(["r"])
    writer.writerow(rZ)'''

fig, ax = plt.subplots()
#ax.scatter(freqs,fitX,label='Eq Signal',color='black')
ax.scatter(freqs,fitHB2,label='HB Signal',color='black')
#ax.scatter(freqs,diff2,label='Difference',color='blue')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

fig, ax = plt.subplots()
#ax.scatter(freqs,simEq,label='Simulated Eq Signal',color='black')
ax.scatter(freqs,cube_subtHB,label='HB Signal',color='black')
ax.plot(freqs,fitHB2,label='HB Fit',color='grey')
ax.plot(freqs,eqComps[0],label='Negative EQ',color='cyan',linestyle='dashed')
ax.plot(freqs,eqComps[1],label='Positive EQ',color='cyan',linestyle='dotted')
ax.plot(freqs,fitComps[0],label='Negative HB',color='red',linestyle='dashed')
ax.plot(freqs,fitComps[1],label='Positive HB',color='red',linestyle='dotted')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()
