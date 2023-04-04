import math
import numpy as np

def lorentzian(freq, f0, w, K, theta, c):
    #w = FWHM, always positive
    #f0 = central frequency
    #K = scaling, neg is flipped
    #c = y intercept
    x = (f0-freq)/(w/2)
    D = K*(x/(1+x**2))
    A = K*(1/(1+x**2))
    R = A*np.cos(theta)-D*np.sin(theta)+c
    # R = A*cos-D*sin
    return R

def lorentz2(freq, f0, gamma, cc, b):
    top = gamma
    bot = (freq-f0)**2 + gamma**2
    out = top/bot
    return cc*out + b

def holeNbump(freq, f0, w, KHB, theta, omegD):
    hole = lorentzian(freq, f0, w, KHB, theta, 0.0)
    bump = lorentzian(freq, 2*omegD - f0, w, -0.5*KHB, theta, 0.0)
    total = hole + bump
    return total

def linear(x, a, b):
    y = a*x + b
    return y

def quadratic(x, a, b, c):
    y = a*x**2 + b*x + c
    return y

def cubic(x, a, b, c, d):
    y = a*x**3 + b*x**2 + c*x + d
    return y

def deut_half(freq, epsi, omegD, omegQ, A, eta, phi):
    #The functional form of the deuteron NMR signal from 'A line-shape analysis for spin-1 NMR signals', Dulya et al.
    #epsi = +/-1
    #phi in radians
    R = (freq - omegD)/(3*omegQ)
    theta = omegQ/omegD
    rhosq = np.sqrt(A**2 + (1 - epsi*R - eta*np.cos(2*phi))**2)
    rho = np.sqrt(rhosq)
    alpha = np.arccos((1 - epsi*R - eta*np.cos(2*phi))/rhosq)
    Y = np.sqrt(3 - eta*np.cos(2*phi))
    left = 2*np.cos(alpha/2)*(np.arctan((Y**2 - rhosq)/(2*Y*rho*np.sin(alpha/2))) + math.pi/2)
    right_num = Y**2 + rhosq + 2*Y*rho*np.cos(alpha/2)
    right_dem = Y**2 + rhosq - 2*Y*rho*np.cos(alpha/2)
    right = np.sin(alpha/2)*np.log(right_num/right_dem)
    unpol = (1/(2*math.pi*rho))*(left + right)  #equation 14 in Dulya
    return unpol

def deut_phiavg(freq, epsi, omegD, omegQ, A, eta, r):
    #averaging the function over phi from 0 to 2pi
    R = (freq - omegD)/(3*omegQ)
    theta = omegQ/omegD
    J = 64
    avgtot = 0.0
    for j in range(J):
        phi = 2*math.pi*j/J
        unavg_num = np.sqrt(3)*deut_half(freq, epsi, omegD, omegQ, A, eta, phi)
        unavg_den = np.sqrt(3 - eta*np.cos(2*phi))
        avgtot += unavg_num/unavg_den
    aver = avgtot/(J+1)
    if epsi == -1:
        aver = aver*((r**(1+3*theta*R) - 1)/(r**(1+theta*R)))
    if epsi == 1:
        aver = aver*((r**2 - r**(1-3*theta*r))/(r**(1-theta*R)))
    abso = aver/omegQ
    return abso

def deut_phiavg_HB(freq, epsi, omegD, omegQ, A, eta, r, f0, w, KHB1, KHB2):
    base = deut_phiavg(freq, epsi, omegD, omegQ, A, eta, r)
    hole = lorentzian(freq, f0, w, KHB1, 0.0, 0.0)
    bump = lorentzian(freq, 2*omegD - f0, w, KHB2, 0.0, 0.0)
    total = base + hole + bump
    '''if (total < 0).any():
        print(str((total < 0).any()))
        total = 0.0'''
    return total

def deut_fit(freq, omegD, omegQ, A, eta, r):
    #both functions added together
    neg = deut_phiavg(freq, -1, omegD, omegQ, A, eta, r)
    pos = deut_phiavg(freq, 1, omegD, omegQ, A, eta, r)
    total = neg + pos    #equation 24 in Dulya
    return total

def deut_fit_HB(freq, omegD, omegQ, A, eta, r, f0, w, Ka, Kb):
    neg = deut_phiavg_HB(freq, -1, omegD, omegQ, A, eta, r, f0, w, Ka, -0.5*Kb)
    pos = deut_phiavg_HB(freq, 1, omegD, omegQ, A, eta, r, f0, w, Kb, -0.5*Ka)
    total = neg + pos
    return total

def deut_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale):
    #carbon and oxygen bonds separately
    A2 = A1*omegQ1/omegQ2
    one = deut_fit(freq, omegD, omegQ1, A1, eta1, r)
    two = deut_fit(freq, omegD, omegQ2, A2, eta2, r)
    total = scale*((1-K)*one + K*two)
    return total

def deut_double_HB(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, f0, w, Ka, Kb):
    A2 = A1*omegQ1/omegQ2
    one = deut_fit_HB(freq, omegD, omegQ1, A1, eta1, r, f0, w, Ka, Kb)
    two = deut_fit_HB(freq, omegD, omegQ2, A2, eta2, r, f0, w, Ka, Kb)
    total = scale*((1-K)*one + K*two)
    return total

def DHB_hardcode(freq, r, f0, w, Ka, Kb):
    omegD = 3.29287000e+01
    omegQ1 = 2.16607624e-02
    omegQ2 = 2.87585358e-02
    A1 = 3.11542336e-02
    eta1 = 6.69345396e-02
    eta2 = 6.73700719e-02
    K = 1.13970002e-01
    scale = 1.10797550e-02
    out = deut_double_HB(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, f0, w, Ka, Kb)
    return out
        
                   
def deut_double_NaiveHB(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, f0, w, KHB, theta):
    base = deut_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale)
    hnb = holeNbump(freq, f0, w, KHB, theta, omegD)
    total = base + hnb
    return total

def deut_disp_half(freq, epsi, omegD, omegQ, A, eta, phi):
    #The functional form of the deuteron NMR signal from 'A line-shape analysis for spin-1 NMR signals', Dulya et al.
    #epsi = +/-1
    #phi in radians
    R = (freq - omegD)/(3*omegQ)
    theta = omegQ/omegD
    rhosq = np.sqrt(A**2 + (1 - epsi*R - eta*np.cos(2*phi))**2)
    rho = np.sqrt(rhosq)
    alpha = np.arccos((1 - epsi*R - eta*np.cos(2*phi))/rhosq)
    Y = np.sqrt(3 - eta*np.cos(2*phi))
    left = 2*np.sin(alpha/2)*(np.arctan((Y**2 - rhosq)/(2*Y*rho*np.sin(alpha/2))) + math.pi/2)
    right_num = Y**2 + rhosq + 2*Y*rho*np.cos(alpha/2)
    right_dem = Y**2 + rhosq - 2*Y*rho*np.cos(alpha/2)
    right = np.cos(alpha/2)*(np.log(right_dem/right_num))
    unpol = 1*(1/(2*math.pi*rho))*(left + right)  #CHANGE
    return unpol

def deut_disp_phiavg(freq, epsi, omegD, omegQ, A, eta, r):
    #averaging the function over phi from 0 to 2pi
    R = (freq - omegD)/(3*omegQ)
    theta = omegQ/omegD
    J = 64
    avgtot = 0.0
    for j in range(J):
        phi = 2*math.pi*j/J
        unavg_num = np.sqrt(3)*deut_disp_half(freq, epsi, omegD, omegQ, A, eta, phi)
        unavg_den = np.sqrt(3 - eta*np.cos(2*phi))
        avgtot += unavg_num/unavg_den
    aver = avgtot/(J+1)
    if epsi == -1:
        aver = -1*aver*((r**(1+3*theta*R) - 1)/(r**(1+theta*R)))
    if epsi == 1:
        aver = aver*((r**2 - r**(1-3*theta*r))/(r**(1-theta*R)))
    abso = aver/omegQ
    return abso

def deut_disp_fit(freq, omegD, omegQ, A, eta, r):
    #both functions added together
    neg = deut_disp_phiavg(freq, -1, omegD, omegQ, A, eta, r)
    pos = deut_disp_phiavg(freq, 1, omegD, omegQ, A, eta, r)
    total = neg + pos    #equation 24 in Dulya
    return total

def deut_disp_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale):
    #carbon and oxygen bonds separately
    A2 = A1*omegQ1/omegQ2
    one = deut_disp_fit(freq, omegD, omegQ1, A1, eta1, r)
    two = deut_disp_fit(freq, omegD, omegQ2, A2, eta2, r)
    total = scale*((1-K)*one + K*two)
    return total

def deut_real_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, phase):
    absor = deut_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale)
    disper = deut_disp_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale)
    real = absor*np.cos(phase) - disper*np.sin(phase)
    return real

def deut_real_double_cube(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, phase, x0, x1, x2, x3):
    recur = deut_real_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, phase)
    cube = cubic(freq, x3, x2, x1, x0)
    return recur + cube

def deut_real_double_lorz(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, phase, f0, gamma, cc, b):
    recur = deut_real_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, phase)
    lorz = lorentz2(freq, f0, gamma, cc, b)
    return recur + lorz

def deut_im_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, phase, offs, slope):
    absor = deut_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale)
    disper = deut_disp_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale)
    imag = absor*np.sin(phase) + disper*np.cos(phase) + offs + slope*freq
    return imag

def DIDhardCode(freq, r, slope, offs):
    omegD = 3.27705415e+01
    omegQ1 = 2.16712844e-02
    omegQ2 = 2.75852632e-02
    A1 = 2.90224589e-02
    eta1 = 1.01142710e-01
    eta2 = 1.37951403e-01
    K = 1.55024437e-01
    scale = 4.84541259e-04
    phase = -3.31813519e-01
    out = deut_im_double(freq,omegD,omegQ1,omegQ2,A1,eta1,eta2,r,K,scale,phase,offs,slope)
    return out

def DRD_fitR(freq, r):
    omegD = 3.27702256e+01
    omegQ1 = 2.16986220e-02
    omegQ2 = 2.59889140e-02
    A1 = 2.97697232e-02
    eta1 = 6.18836853e-02
    eta2 = 1.98712673e-01
    K = 1.94099997e-01
    scale = 5.12571386e-04
    phase = 1.90848493e-03
    out = deut_real_double(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, phase)
    return out

def DRD_fitR_cube(freq, r, x0, x1, x2, x3):
    omegD=3.30048700e+01
    omegQ1=2.16427904e-02
    omegQ2=2.70776767e-02
    A1=2.61510356e-02
    eta1=6.62185536e-02
    eta2=1.46435235e-01
    K=9.75343811e-02
    scale=2.69027249e-03
    phase=3.35679145e-02
    recur = deut_real_double_cube(freq, omegD, omegQ1, omegQ2, A1, eta1, eta2, r, K, scale, phase, x0, x1, x2, x3)
    return recur

def ratio_method(r):
    pol = (r*r - 1)/(r + r*r + 1)
    return pol

def ratio_tensor(r):
    pol = (r*r - 2*r + 1)/(r + r*r + 1)
    return pol

def riemann_step(ydata, step):
    area = 0.0
    Range = len(ydata)
    for i in range(Range):
        area += step*ydata[i]
    return area

def riemann_sum(xdata, ydata):
    area = 0.0
    Range = len(xdata)
    minim = xdata[0]
    maxim = xdata[Range-1]
    step = (maxim - minim)/Range
    for i in range(Range):
        area += step*ydata[i]
    return area

def freq2R(freq):
    R = 16.32516*freq - 537.58733
    return R
