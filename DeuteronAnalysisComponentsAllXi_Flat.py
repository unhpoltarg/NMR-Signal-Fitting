import math
import numpy as np
import matplotlib.pyplot as plt
import csv
from datetime import datetime
from scipy.optimize import curve_fit
from datetime import datetime
import deuts
import time

'''
Your data file should match ExampleCSV.csv:
Line 1: Central Frequency | Frequency Span | NMR Width (number of data points)
Line 2: Baseline data
Line 3: Signal data
***
'''
#YOU MUST HAVE THE deuts.py FILE IN THE SAME FOLDER

def DRDC_HC(freq,omegD,r,scale,phase,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase,x3,x2,x1,x0)
    return out

newRaw = [] #Your data
freqs = []

centfreq = 0.0
freqspan = 0.0
nmrwidth = 0

filename = 'ExampleFlatCSV.csv' #Change to your csv's name - should be in same folder as this file
file = open(filename, 'r')
reader = csv.reader(file)

rownum = 1
for row in reader:
    if rownum == 1:
        nmrwidth = int(row[0])
    if rownum == 2:
        for i in range(nmrwidth):
            freqs.append(float(row[0+i]))
    if rownum == 3:
        for i in range(nmrwidth):
            newRaw.append(float(row[0+i]))
    rownum += 1

file.close()
centfreq = freqs[int(nmrwidth/2)]

print("Cent freq = " + str(centfreq))
print(newRaw[0])

#if signal is too wide, comment this in
'''tempF = []
tempD = []
for i in range(nmrwidth):
    if i > 60 and i < nmrwidth - 60:
        tempF.append(freqs[i])
        tempD.append(newRaw[i])
freqs = tempF
newRaw = tempD
nmrwidth = len(freqs)'''
        

#This will produce a plot of your flattened signal. If the signal
# looks wrong, go back and check your csv to make sure it adheres to requirements
fig, ax = plt.subplots()
ax.plot(freqs,newRaw,label='Signal',color='black')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

#This fits the flattened signal. The parameters will print out in IDLE.
#If any parameters are "railing" (going to X.99999 or X.00001), adjust the minima or maxima.
pFit, cFit = curve_fit(deuts.deut_real_double_Xi, freqs, newRaw, maxfev=10000,
                       p0=[centfreq,0.021711,0.026,0.026,0.06,0.2,1.5,0.08,0.004019,0],
                       bounds=([centfreq-0.1,0.0212,0.025,0.0001,0.0001,0.15,0.1,0.03,0,-2.5],
                               [centfreq+0.1,0.022,0.028,0.1,0.12,0.25,3,0.25,0.02,2.5]))
errors = np.sqrt(np.diag(cFit))
print(np.sqrt(abs(cFit[2][4])))
print(np.sqrt(cFit[2][2]))
#use this if already flat
'''pFit, cFit = curve_fit(deuts.deut_real_double_Xi, freqs, newRaw, maxfev=100000,
                           bounds=([centfreq-0.1,0.0212,0.027,0.01,0.01,0.01,1,0.01,0,-0.5],
                               [centfreq+0.1,0.0218,0.031,0.04,0.1,0.2,3,0.3,0.05,0.5]))'''


omegD = pFit[0] # Omega_D - the central frequency of the signal
print("Larmor frequency (omegaD) = " + format(omegD, '.6f') + " ± " + format(errors[0], '.6f'))
# Omega_Q - the quadrupole frequency: 3*omegaQ = distance from omegaD to the peak
omegQ1 = pFit[1] # Omega_Q1 - carbon quadrupole frequency
print("Carbon quadrupole frequency (omegaQ1) = " + format(omegQ1, '.6f') + " ± " + format(errors[1], '.6f'))
omegQ2 = pFit[2] # Omega_Q2 - oxygen quadrupole frequency
print("Oxygen quadrupole frequency (omegaQ2) = " + format(omegQ2, '.6f') + " ± " + format(errors[2], '.6f'))
# A - the Lorentzian width of the signal
A1 = pFit[3] # A1 - carbon Lorentzian width (oxygen's is pegged to it: omegaQ1*A1 = omegaQ2*A2)
print("Carbon Lorentzian width (A1) = " + format(A1, '.6f') + " ± " + format(errors[3], '.6f'))
A2 = (A1*omegQ1)/omegQ2
# eta - asymmetry parameter
eta1 = pFit[4] # eta1 - carbon
print("Carbon electric asymmetry (eta1) = " + format(eta1, '.6f') + " ± " + format(errors[4], '.6f'))
eta2 = pFit[5] # eta2 - oxygen
print("Oxygen electric asymmetry (eta2) = " + format(eta2, '.6f') + " ± " + format(errors[5], '.6f'))
r = pFit[6] # r - "ratio" number that describes Pz and Pzz
print("Signal asymmetry (r) = " + format(r, '.6f') + " ± " + format(errors[6], '.6f'))
K = pFit[7] # K - the percentage of deuteron bonds that are with oxygen
print("Oxygen proportion (K) = " + format(K, '.6f') + " ± " + format(errors[7], '.6f'))
scale = pFit[8] # scaling factor - turns the proportionality of chi'' into a function
print("Constant scaling (c) = " + format(scale, '.6f') + " ± " + format(errors[8], '.6f'))
print(scale)
phase = pFit[9] # phase angle - if the signal is complex
print("Phase angle (theta) = " + format((phase*180)/np.pi, '.6f') + " ± " + format((errors[9]*180)/np.pi, '.6f') + " degrees")


#designate the fit
fit = []
data_std = 0.0
for i in range(nmrwidth):
    fit.append(deuts.deut_real_double_Xi(freqs[i],omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase))
    data_std += abs(newRaw[i] - fit[i])
data_std = data_std/nmrwidth
print("Standard deviation of data = " + str(data_std))

Perr = 100*deuts.ratioPderiv(r)*errors[6]
Qerr = 100*deuts.ratioQderiv(r)*errors[6]

#Print in IDLE Pz and Pzz in decimal form.
Pz = deuts.ratio_method(r)
Pzz = deuts.ratio_tensor(r)
print("Fit P = " + format(100*Pz, '.6f') + " ± " + format(Perr, '.6f') + "%")
print("Fit Q = " + format(100*Pzz, '.6f') + " ± " + format(Qerr, '.6f') + "%")
STD = 0.0009863614666752998
print("Reduced chi-squared = " + str(deuts.chisq_red(newRaw,fit,STD)))
print("Reduced chi-squared from data = " + str(deuts.chisq_red(newRaw,fit,data_std)))
area = deuts.riemann_sum(freqs,fit)
CC = Pz/area
print("CC = " + str(CC))

AreaPz = CC*area
print("Area P = " + str(100*AreaPz)+"%")

#component spin-flip curves
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
    out = deuts.Xi(r)*scale*k*(deuts.deut_phiavg(freq,epsi,omegD,omegQ,A,eta,r)*np.cos(phase) + 
            deuts.deut_disp_phiavg(freq,epsi,omegD,omegQ,A,eta,r)*np.sin(phase))
    return out

#Get component curves of fit
fitComps = [[] for i in range(4)]
for i in range(nmrwidth):
    fitComps[0].append(NormComp(freqs[i],'C',-1))
    fitComps[1].append(NormComp(freqs[i],'C',1))
    fitComps[2].append(NormComp(freqs[i],'O',-1))
    fitComps[3].append(NormComp(freqs[i],'O',1))

Tarea = deuts.riemann_sum(freqs,fitComps[1]) - deuts.riemann_sum(freqs,fitComps[0]) + deuts.riemann_sum(freqs,fitComps[3]) - deuts.riemann_sum(freqs,fitComps[2])
Parea = deuts.riemann_sum(freqs,fitComps[1]) + deuts.riemann_sum(freqs,fitComps[3]) + deuts.riemann_sum(freqs,fitComps[0]) + deuts.riemann_sum(freqs,fitComps[2])

AreaP = CC*Parea
AreaPzz = CC*Tarea
print("Area Sum P = " + str(100*AreaP)+"%")
print("Area Q = " + str(100*AreaPzz)+"%")
resid = []
onlyOx = []
fitOx = []
for i in range(nmrwidth):
    resid.append(fit[i]-newRaw[i])
    onlyOx.append(newRaw[i] - fitComps[0][i] - fitComps[1][i])
    fitOx.append(fitComps[2][i] + fitComps[3][i])

#This shows a plot of your fit compared to the data.
fig, ax = plt.subplots()
ax.scatter(freqs,newRaw,label='Data',color='black')
ax.plot(freqs,fit,label='Fit',color='gray')
ax.plot(freqs,fitComps[0],label='C-D Negative',color='red',linestyle='dashed')
ax.plot(freqs,fitComps[1],label='C-D Positive',color='red',linestyle='dotted')
ax.plot(freqs,fitComps[2],label='O-D Negative',color='blue',linestyle='dashed')
ax.plot(freqs,fitComps[3],label='O-D Positive',color='blue',linestyle='dotted')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

fig, ax = plt.subplots()
ax.scatter(freqs,resid,label='Data',color='black')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Residual (a.u.)")
plt.legend()
plt.show() 
