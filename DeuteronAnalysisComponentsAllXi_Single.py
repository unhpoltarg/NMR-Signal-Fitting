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

#expected values for propanediol - these will need to be adjusted for other materials
#omegQ1 = 2.17139774e-02
#omegQ2 = 2.88692494e-02
#A1 = 3.06372826e-02
#eta1 = 4.59723046e-02
#eta2 = 6.12644049e-02
#K = 1.15523343e-01

def DRDC_HC(freq,omegD,r,scale,phase,x3,x2,x1,x0):
    out = deuts.deut_real_double_cube_Xi(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase,x3,x2,x1,x0)
    return out

rawdata = [] #Your data
baseline = [] #your baseline
bases = []

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
            baseline.append(-1*float(row[0+i]))
    if rownum == 3:
        for i in range(nmrwidth):
            rawdata.append(-1*float(row[0+i]))
    '''if rownum >= 5:
        temp = []
        for i in range(nmrwidth):
            temp.append(float(row[0+i]))
        bases.append(temp)'''
    rownum += 1

file.close()

'''avgbase = []
for i in range(nmrwidth):
    temp = 0.0
    for j in range(len(bases)):
        temp += bases[j][i]
    temp = temp/len(bases)
    avgbase.append(temp)'''

print("Cent freq = " + str(centfreq))

freqs = [] #This propogates the frequencies associated with your data points

for i in range(nmrwidth):
    freqs.append(centfreq + freqspan*(i/(nmrwidth*1.0) - 0.5))

newRaw = [] #Subtracting the baseline from the signal

for i in range(nmrwidth):
    newRaw.append(rawdata[i] - baseline[i])#baseline[i])

freqs = freqs[80:nmrwidth - 80]
newRaw = newRaw[80:nmrwidth - 80]

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
pFit, cFit = curve_fit(deuts.deut_real_cube_Xi, freqs, newRaw, maxfev=100000,
                           bounds=([centfreq-0.05,0.0214,0.01,0.001,1,0,-0.5,-100,-1000,-1000,-100000],
                               [centfreq+0.05,0.03,0.5,0.5,3,0.1,0.5,100,1000,1000,100000]))

errors = np.sqrt(np.diag(cFit))
omegD = pFit[0] # Omega_D - the central frequency of the signal
print("Larmor frequency (omegaD) = " + format(omegD, '.6f') + " ± " + format(errors[0], '.6f'))
# Omega_Q - the quadrupole frequency: 3*omegaQ = distance from omegaD to the peak
omegQ1 = pFit[1] # Omega_Q1 - carbon quadrupole frequency
print("Quadrupole frequency (omegaQ) = " + format(omegQ1, '.6f') + " ± " + format(errors[1], '.6f'))
# A - the Lorentzian width of the signal
A1 = pFit[2] # A1 - carbon Lorentzian width (oxygen's is pegged to it: omegaQ1*A1 = omegaQ2*A2)
print("Lorentzian width (A) = " + format(A1, '.6f') + " ± " + format(errors[2], '.6f'))
# eta - asymmetry parameter
eta1 = pFit[3] # eta1 - carbon
print("Electric asymmetry (eta) = " + format(eta1, '.6f') + " ± " + format(errors[3], '.6f'))
r = pFit[4] # r - "ratio" number that describes Pz and Pzz
print("Signal asymmetry (r) = " + format(r, '.6f') + " ± " + format(errors[4], '.6f'))
scale = pFit[5] # scaling factor - turns the proportionality of chi'' into a function
print("Constant scaling (c) = " + format(scale, '.6f') + " ± " + format(errors[5], '.6f'))
phase = pFit[6] # phase angle - if the signal is complex
print("Phase angle (theta) = " + format((phase*180)/np.pi, '.6f') + " ± " + format((errors[6]*180)/np.pi, '.6f') + " degrees")

# cubic subtraction
x3 = pFit[7]
print(x3)
x2 = pFit[8]
print(x2)
x1 = pFit[9]
print(x1)
x0 = pFit[10]
print(x0)
#designate the fit
coobic = []
cube_subt = []
fit = []
for i in range(len(freqs)):
    coobic.append(deuts.cubic(freqs[i],x3,x2,x1,x0))
    cube_subt.append(newRaw[i] - coobic[i])
    fit.append(deuts.deut_real_Xi(freqs[i],omegD,omegQ1,A1,eta1,r,scale,phase))

fig, ax = plt.subplots()
ax.scatter(freqs,newRaw,label='Signal',color='black')
ax.scatter(freqs,coobic,label='Cubic',color='blue')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show()

Perr = 100*deuts.ratioPderiv(r)*errors[4]
Qerr = 100*deuts.ratioQderiv(r)*errors[4]

#Print in IDLE Pz and Pzz in decimal form.
Pz = deuts.ratio_method(r)
Pzz = deuts.ratio_tensor(r)
print("Fit P = " + format(100*Pz, '.6f') + " ± " + format(Perr, '.6f') + "%")
print("Fit Q = " + format(100*Pzz, '.6f') + " ± " + format(Qerr, '.6f') + "%")
STD = 0.0009863614666752998
print("Reduced chi-squared = " + str(deuts.chisq_red(cube_subt,fit,STD)))

#Get component curves of fit
fitComps = [[] for i in range(4)]
for i in range(len(freqs)):
    fitComps[0].append(deuts.Xi(r)*scale*(deuts.deut_phiavg(freqs[i],-1,omegD,omegQ1,A1,eta1,r)*np.cos(phase) + 
        deuts.deut_disp_phiavg(freqs[i],-1,omegD,omegQ1,A1,eta1,r)*np.sin(phase))) #Negative carbon
    fitComps[1].append(deuts.Xi(r)*scale*(deuts.deut_phiavg(freqs[i],1,omegD,omegQ1,A1,eta1,r)*np.cos(phase) + 
        deuts.deut_disp_phiavg(freqs[i],1,omegD,omegQ1,A1,eta1,r)*np.sin(phase))) #Positive carbon

#This shows a plot of your fit compared to the data.
fig, ax = plt.subplots()
ax.plot(freqs,cube_subt,label='Data',color='black')
ax.plot(freqs,fit,label='Fit',color='grey')
ax.plot(freqs,fitComps[0],label='Negative',color='red',linestyle='dashed')
ax.plot(freqs,fitComps[1],label='Positive',color='red',linestyle='dotted')
plt.xlabel("Frequency (MHz)")
plt.ylabel("Signal (a.u.)")
plt.legend()
plt.show() 
