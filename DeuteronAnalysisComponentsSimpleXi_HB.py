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
omegQ1 = 0.021729
omegQ2 = 0.025147
A1 = 0.021409
A2 = (A1*omegQ1)/omegQ2
eta1 = 0.033559
eta2 = 0.230762
K = 0.082220
scale = 0.000326

#expected values for TEMPO butanol
'''omegQ1 = 2.17004538e-02
omegQ2 = 2.89859414e-02
A1 = 2.75141728e-02
eta1 = 6.72464295e-02
eta2 = 6.57934478e-02
K = 4.28479420e-02'''

#expected values for irradiated butanol
'''omegQ1 = 0.021689
omegQ2 = 0.026200
A1 = 0.025325
A2 = (A1*omegQ1)/omegQ2
eta1 = 0.036339
eta2 = 0.174368
K = 0.070193
scale = 0.000064'''
#K = 0.136 according to Dulya, seems wrong
#scale = 9.66383779e-04

'''def DRDC_HC(freq,omegD,r,scale,phase,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase,x3,x2,x1,x0)
    return out

def DRDC_Xi_HC(freq,omegD,r,scale,phase,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase,x3,x2,x1,x0)
    return out'''

def DRDC_XiO_HC(freq,omegD,r,phaseX,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phaseX,x3,x2,x1,x0)
    return out

def DRDC_XiO_scHC(freq,omegD,r,scaleX,phaseX,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scaleX,phaseX,x3,x2,x1,x0)
    return out

def DRD_XiO_scHC(freq,omegD,r,scaleX,phaseX):
    out = deuts.deut_real_double_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scaleX,phaseX)
    return out

def DRD_XiO_HC_NC(freq,omegD,r,phaseX):
    out = deuts.deut_real_double_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phaseX)
    return out

print(str(deuts.ratio_method(2.5)))

rawdata = [] #Your equilibrium data
HBdata = [] #Your holeburn data
baseline = [] #your baseline
holeline = []

centfreq = 0.0
freqspan = 0.0
nmrwidth = 0

filename = 'ExampleCSV_HB1.csv' #Change to your csv's name - should be in same folder as this file
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
    if rownum == 4:
        for i in range(nmrwidth):
            HBdata.append(float(row[0+i]))
    rownum += 1

file.close()

print("Cent freq = " + str(centfreq))

freqs = [] #This propogates the frequencies associated with your data points

for i in range(nmrwidth):
    freqs.append(centfreq + freqspan*(i/(nmrwidth*1.0) - 0.5))

newRaw = [] #Subtracting the baseline from the signals
newHB = []

for i in range(nmrwidth):
    newRaw.append(rawdata[i] - baseline[i])
    newHB.append(HBdata[i] - baseline[i])# - holeline[i])#

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
                           bounds=([centfreq-0.5,1,-3,-1000,-1000,-10000,-100000],
                               [centfreq+0.5,3,3,1000,1000,10000,100000]))

print(pFit2)
omegD = pFit2[0]
rX = pFit2[1]
#OVERRIDE
#rX = 1.08
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
    fitX.append(deuts.deut_real_double_Xi(freqs[i],omegD,omegQ1,omegQ2,A1,eta1,eta2,rX,K,scale,phase))

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
    if i < 25 or i > nmrwidth - 25:
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
#ax.plot(freqs,fitX,label='Fit',color='red')
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
    outt = deuts.deut_real_double_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,arr,K,scale,phase)
    return outt

fitHB = []
fitHB_IM = []
fitComps = [[] for i in range(4)]
fitCompsIm = [[] for i in range(4)]
eqComps = [[] for i in range(4)]
eqCompsIm = [[] for i in range(4)]

def Bump(freq,req,arr):
            outt = (deuts.deut_real_double_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,req,K,scale,phase) +
                    0.5*deuts.deut_real_double_Xi_rev(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,req,K,scale,phase) -
                    0.5*deuts.deut_real_double_Xi_rev(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,arr,K,scale,phase))                                                                                                        
            return outt

def NormComp(freq,arr,CO,phi,epsi):
    if CO == 'C':
        k = 1 - K
        omegQ = omegQ1
        A = A1
        eta = eta1
    if CO == 'O':
        k = K
        omegQ = omegQ2
        A = A2
        eta = eta2
    out = deuts.Xi(arr)*scale*k*(deuts.deut_phiavg(freq,epsi,omegD,omegQ,A,eta,arr)*np.cos(phi) + 
            deuts.deut_disp_phiavg(freq,epsi,omegD,omegQ,A,eta,arr)*np.sin(phi))
    return out

def NormCompIm(freq,arr,CO,epsi):
    if CO == 'C':
        k = 1 - K
        omegQ = omegQ1
        A = A1
        eta = eta1
    if CO == 'O':
        k = K
        omegQ = omegQ2
        A = A2
        eta = eta2
    out = deuts.Xi(arr)*scale*k*(deuts.deut_phiavg(freq,epsi,omegD,omegQ,A,eta,arr)*np.cos(math.pi/2.0) + 
            deuts.deut_disp_phiavg(freq,epsi,omegD,omegQ,A,eta,arr)*np.sin(math.pi/2.0))
    return out

def BumpComp(freq,req,rstc,CO,phi,epsi): #CO = C or O
    if CO == 'C':
        k = 1 - K
        omegQ = omegQ1
        A = A1
        eta = eta1
    if CO == 'O':
        k = K
        omegQ = omegQ2
        A = A2
        eta = eta2
    equi = NormComp(freq,req,CO,phi,epsi)
    fake = 0.5*deuts.Xi(req)*scale*k*(deuts.deut_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,req)*np.cos(phi) + 
                                deuts.deut_disp_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,req)*np.sin(phi))
    holeSub = 0.5*deuts.Xi(rstc)*scale*k*(deuts.deut_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,rstc)*np.cos(phi) + 
                                deuts.deut_disp_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,rstc)*np.sin(phi))
    out = equi + fake - holeSub
    return out

def BumpCompIm(freq,req,rstc,CO,epsi): #CO = C or O
    if CO == 'C':
        k = 1 - K
        omegQ = omegQ1
        A = A1
        eta = eta1
    if CO == 'O':
        k = K
        omegQ = omegQ2
        A = A2
        eta = eta2
    equi = NormCompIm(freq,req,CO,epsi)
    fake = 0.5*deuts.Xi(req)*scale*k*(deuts.deut_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,req)*np.cos(math.pi/2.0) + 
                                deuts.deut_disp_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,req)*np.sin(math.pi/2.0))
    holeSub = 0.5*deuts.Xi(rstc)*scale*k*(deuts.deut_phiavg_rev(freq,epsi,omegD,omegQ,A,eta,rstc)*np.cos(math.pi/2.0) + 
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

Pbins = []
Qbins = []
bin87 = []
for i in range(nmrwidth):
    Pbins.append(100*deuts.ratio_method(rZ[i]))
    Qbins.append(100*deuts.ratio_tensor(rZ[i]))
    bin87.append(DRD_XiO_HC_NC(freqs[i],omegD,rZ[87],phase))

fig, ax = plt.subplots()
ax.scatter(freqs,cube_subtHB,label='Data',color='black')
ax.plot(freqs,bin87,label="Fit",color='blue')
ax.plot(freqs[87],cube_subtHB[87],'ro',label="Point")
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

fig, ax = plt.subplots()
ax.scatter(freqs,rZ,label='r',color='green')
plt.xlabel("Frequency (MHz)")
plt.ylabel("P")
plt.legend()
plt.show()

rBump.append(1.0)
fitHB = []
for i in range(len(freqs)):
        fitHB.append(deuts.deut_real_double_Xi(freqs[i],omegD,omegQ1,omegQ2,A1,eta1,eta2,rZ[i],K,scale,phase))
        #fitHB_IM.append(deuts.deut_real_double_Xi(freqs[i],omegD,omegQ1,omegQ2,A1,eta1,eta2,rZ[i],K,scaleX,math.pi/2.0))

for i in range(nmrwidth):
    eqComps[0].append(NormComp(freqs[i],rX,'C',phase,-1))
    eqComps[1].append(NormComp(freqs[i],rX,'C',phase,1))
    eqComps[2].append(NormComp(freqs[i],rX,'O',phase,-1))
    eqComps[3].append(NormComp(freqs[i],rX,'O',phase,1))
    if rZ[i] > rX:
        fitComps[0].append(BumpComp(freqs[i],rX,rBump[i],'C',phase,-1))
        fitComps[1].append(BumpComp(freqs[i],rX,rBump[i],'C',phase,1))
        fitComps[2].append(BumpComp(freqs[i],rX,rBump[i],'O',phase,-1))
        fitComps[3].append(BumpComp(freqs[i],rX,rBump[i],'O',phase,1))
            
    else:
        fitComps[0].append(NormComp(freqs[i],rZ[i],'C',phase,-1))
        fitComps[1].append(NormComp(freqs[i],rZ[i],'C',phase,1))
        fitComps[2].append(NormComp(freqs[i],rZ[i],'O',phase,-1))
        fitComps[3].append(NormComp(freqs[i],rZ[i],'O',phase,1))

#simEq = []
fitHB2 = []
for i in range(len(rZ)):
    fitHB2.append(deuts.deut_real_double_Xi(freqs[i],omegD,omegQ1,omegQ2,A1,eta1,eta2,rZ[i],K,scale,phase))
    #simEq.append(deuts.deut_real_double_Xi(freqs[i],omegDX,omegQ1,omegQ2,A1,eta1,eta2,rZ[226],K,scaleX,phaseHB))

#Print in IDLE Pz and Pzz in decimal form.
STD = 0.0009863614666752998

PzX = deuts.ratio_method(rX)
PzzX = deuts.ratio_tensor(rX)
print("P = " + str(100*PzX)+"%")
print("Q = " + str(100*PzzX)+"%")
print("Phase angle = " + str((180*phase)/np.pi) + " degrees")
print("Reduced chi-squared = " + str(deuts.chisq_red(cube_subtX,fitX,STD)))
CCeq = PzX/deuts.riemann_sum(freqs,cube_subtX)

eqSTD = 0.0
for i in range(nmrwidth):
    eqSTD += abs(fitX[i] - cube_subtX[i])
eqSTD = (eqSTD/nmrwidth)*CCeq*100
print("Error = " + str(eqSTD)+"%")
print("CC = " + str(CCeq))
EminC = deuts.riemann_sum(freqs,eqComps[0])
EplusC = deuts.riemann_sum(freqs,eqComps[1])
EminO = deuts.riemann_sum(freqs,eqComps[2])
EplusO = deuts.riemann_sum(freqs,eqComps[3])
IminC = deuts.riemann_sum(freqs,fitComps[0])
IplusC = deuts.riemann_sum(freqs,fitComps[1])
IminO = deuts.riemann_sum(freqs,fitComps[2])
IplusO = deuts.riemann_sum(freqs,fitComps[3])
EQA_Q = CCeq*(EplusC + EplusO - EminC - EminO)
HB_P = CCeq*(IplusC + IplusO + IminC + IminO)
HB_Q = CCeq*(IplusC + IplusO - IminC - IminO)
larm = int(nmrwidth/2)
P_larm = deuts.ratio_method(rZ[larm])
Q_larm = deuts.ratio_tensor(rZ[larm])
cyanTot = eqComps[0][larm] + eqComps[1][larm] + eqComps[2][larm] + eqComps[3][larm]
greyTot = fitHB2[larm]
print("P larmor = " + str(100*P_larm))
print("Q larmor = " + str(100*Q_larm))
print("EQ Larmor = " + str(cyanTot))
print("HB Larmor = " + str(greyTot))
print("EQ I plus = " + str(CCeq*(EplusC + EplusO)))
print("EQ I minus = " + str(CCeq*(EminC + EminO)))
print("HB I plus = " + str(CCeq*(IplusC + IplusO)))
print("HB I minus = " + str(CCeq*(IminC + IminO)))
    
print("Area Q = "  +str(100*EQA_Q)+"%")
print("Hole burnt P = " +str(100*HB_P)+"%")
print("Hole burnt Q = " +str(100*HB_Q)+"%")

diff2 = []
for i in range(nmrwidth):
    diff2.append(fitHB2[i] - fitX[i])

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
#ax.plot(freqs,cube_subtX,label='EQ Fit',color='cyan')
#ax.plot(freqs,eqComps[0],label='C-D Negative EQ',color='cyan',linestyle='dashed')
#ax.plot(freqs,eqComps[1],label='C-D Positive EQ',color='cyan',linestyle='dotted')
ax.plot(freqs,fitComps[0],label='C-D Negative',color='red',linestyle='dashed')
ax.plot(freqs,fitComps[1],label='C-D Positive',color='red',linestyle='dotted')
ax.plot(freqs,fitComps[2],label='O-D Negative',color='blue',linestyle='dashed')
ax.plot(freqs,fitComps[3],label='O-D Positive',color='blue',linestyle='dotted')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()
