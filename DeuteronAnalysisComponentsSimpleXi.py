import math
import numpy as np
import matplotlib.pyplot as plt
import csv
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
omegQ1 = 0.021572
omegQ2 = 0.027101
A1 = 0.031919
A2 = (A1*omegQ1)/omegQ2
eta1 = 0.062547
eta2 = 0.147506
K = 0.169736
scale = 0.005071
#phase = 0.951243*(np.pi/180)

#expected values for TEMPO butanol
'''omegQ1 = 2.17004538e-02
omegQ2 = 2.89859414e-02
A1 = 2.75141728e-02
eta1 = 6.72464295e-02
eta2 = 6.57934478e-02
K = 4.28479420e-02'''

#expected values for irradiated butanol
'''omegQ1 = 0.021725
omegQ2 = 0.026224
A1 = 0.024600
A2 = (A1*omegQ1)/omegQ2
eta1 = 0.035014
eta2 = 0.179386
K = 0.070492
scale = 0.000843'''
#K = 0.136 according to Dulya, seems wrong

def DRDC_HC(freq,omegD,r,scale,phase,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase,x3,x2,x1,x0)
    return out

def DRDC_Xi_HC(freq,omegD,r,scale,phase,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase,x3,x2,x1,x0)
    return out

def DRDC_XiO_HC(freq,omegD,r,phase,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase,x3,x2,x1,x0)
    return out

rawdata = [] #Your data
baseline = [] #your baseline

centfreq = 0.0
freqspan = 0.0
nmrwidth = 0

filename = 'ExampleCSV1.csv' #Change to your csv's name - should be in same folder as this file
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

print("Cent freq = " + str(centfreq))

freqs = [] #This propogates the frequencies associated with your data points

for i in range(nmrwidth):
    freqs.append(centfreq + freqspan*(i/(nmrwidth*1.0) - 0.5))

newRaw = [] #Subtracting the baseline from the signal

for i in range(nmrwidth):
    newRaw.append(rawdata[i] - baseline[i]) #

#This will produce a plot of your flattened signal. If the signal
# looks wrong, go back and check your csv to make sure it adheres to requirements.
fig, ax = plt.subplots()
ax.scatter(freqs,newRaw,label='Signal',color='black')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

#This fits the flattened signal. The parameters will print out in IDLE.
#If any parameters are "railing" (going to X.99999 or X.00001), adjust the minima or maxima.
pFit2, cFit2 = curve_fit(DRDC_XiO_HC, freqs, newRaw, maxfev=100000,
                           bounds=([centfreq-0.1,0.1,-2,-1000,-1000,-10000,-100000],
                               [centfreq+0.1,3,2,1000,1000,10000,100000]))

print(pFit2)
omegDX = pFit2[0]
rX = pFit2[1]
phaseX = pFit2[2]

x3x = pFit2[3]
x2x = pFit2[4]
x1x = pFit2[5]
x0x = pFit2[6]

cube_subtX = []
fitX = []
for i in range(nmrwidth):
    cube_subtX.append(newRaw[i] - deuts.cubic(freqs[i],x3x,x2x,x1x,x0x))
    fitX.append(deuts.deut_real_double_Xi(freqs[i],omegDX,omegQ1,omegQ2,A1,eta1,eta2,rX,K,scale,phaseX))

errors = np.sqrt(np.diag(cFit2))
Perr = 100*deuts.ratioPderiv(rX)*errors[1]
Qerr = 100*deuts.ratioQderiv(rX)*errors[1]

#Print in IDLE Pz and Pzz in decimal form.
STD = 0.0009863614666752998
PzX = deuts.ratio_method(rX)
PzzX = deuts.ratio_tensor(rX)
print("Fit P = " + format(100*PzX, '.6f') + " ± " + format(Perr, '.6f') + "%")
print("Fit Q = " + format(100*PzzX, '.6f') + " ± " + format(Qerr, '.6f') + "%")
print("Phase angle (theta) = " + format((phaseX*180)/np.pi, '.6f') + " ± " + format((errors[2]*180)/np.pi, '.8f') + " degrees")
print("Reduced chi-squared = " + str(deuts.chisq_red(cube_subtX,fitX,STD)))
CCX = PzX/deuts.riemann_sum(freqs,fitX)
print("CC = " + str(CCX))

def NormComp(freq,CO,epsi):
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
    out = deuts.Xi(rX)*scale*k*(deuts.deut_phiavg(freq,epsi,omegDX,omegQ,A,eta,rX)*np.cos(phaseX) + 
            deuts.deut_disp_phiavg(freq,epsi,omegDX,omegQ,A,eta,rX)*np.sin(phaseX))
    return out

#Get component curves of fit
fitComps = [[] for i in range(4)]
for i in range(nmrwidth):
    fitComps[0].append(NormComp(freqs[i],'C',-1))
    fitComps[1].append(NormComp(freqs[i],'C',1))
    fitComps[2].append(NormComp(freqs[i],'O',-1))
    fitComps[3].append(NormComp(freqs[i],'O',1))
    
#This shows a plot of your fit compared to the data.
fig, ax = plt.subplots()
ax.scatter(freqs,cube_subtX,label='Data',color='black')
ax.plot(freqs,fitX,label='Fit',color='gray')
#ax.plot(freqs,fitX,label='Xi Fit',color='pink',linestyle='dashed')
ax.plot(freqs,fitComps[0],label='C-D Negative',color='red',linestyle='dashed')
ax.plot(freqs,fitComps[1],label='C-D Positive',color='red',linestyle='dotted')
ax.plot(freqs,fitComps[2],label='O-D Negative',color='blue',linestyle=(0, (3, 1, 1, 1)))
ax.plot(freqs,fitComps[3],label='O-D Positive',color='blue',linestyle=(0, (3, 5, 1, 5, 1, 5)))
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show() 
