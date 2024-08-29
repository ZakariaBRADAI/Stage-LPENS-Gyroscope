from GyroSolver import *
from GyroPhysicalQTE import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def AngleRef(t, pulsation, angle0):
    return pulsation * t + angle0


def Sinusoide(t, pulsation, amplitude, phase0):
    return amplitude * np.cos(pulsation * t + phase0)

def Sinusoide_Standard(t, pulsation, phase0):
    return np.cos(pulsation * t + phase0)

def Get_Mvt_Freq(t, angle, omega_hint, angle0_hint = 0, full_output='False'):
    tab, pcov = curve_fit(AngleRef, t, angle, p0=[omega_hint, angle0_hint])
    match full_output:
        case 'False':
            return tab[0], np.sqrt(pcov[0, 0])
        case 'True':
            return tab[0], np.sqrt(pcov[0, 0]), tab[1], np.sqrt(pcov[1, 1])


def Get_Larmor_Freq(t, angle, omega_hint, angle0_hint = 0, full_output='False'):
    return Get_Mvt_Freq(t, angle, omega_hint, angle0_hint, full_output)


def Get_Rabi_Freq(t, theta, params_hint, full_output='False'):
    '''Genius use cos function.'''
    signal = np.cos(theta)
    tab, pcov = curve_fit(Sinusoide, t, signal, p0=params_hint)
    match full_output:
        case 'False':
            return tab[0], np.sqrt(pcov[0, 0])
        case 'True':
            return tab, np.sqrt(np.diag(pcov))


def Rabi_Pulsation_Amplitude_Dependency(tab_x0, tspan, CI, params, forcing='X', precision='LOW'):
    
    list_Omega_R = []
    
    for x0 in tab_x0:
        params[-3] = x0
        t, theta, _, _, _, _, _ = Gyro_Solver(tspan, CI, params, forcing=forcing, precision=precision)
        omega_R_hint = Compute_Carac_Pulsations(params, CI)[1]
        params_hint = [omega_R_hint, 1, CI[0]]
        
        list_Omega_R.append( Get_Rabi_Freq(t, theta, params_hint, full_output='False')[0] )
    list_Omega_R = np.array(list_Omega_R)
    return list_Omega_R
        


def Verif_momentum_varphi_forced(t, theta, p_phi, omega_R, params, forcing='X'):
    _, m, h, _, _, x0, Phi, omega = params
    match forcing:
        case 'X':
            p_phi_th = -0.5 * m * h * x0 * (omega**2) * np.cos(Phi) * np.cos(theta) / omega_R
            label_th = r'$ f(t) = \dfrac{mhx_0 \omega_L^2}{2 \Omega_R}  \cos(\theta) $'
        case 'XY':
            p_phi_th = m * h * x0 * (omega**2) * np.cos(Phi) * np.cos(theta) / omega_R
            label_th = r'$ f(t) = \dfrac{mhx_0 \omega_L^2}{\Omega_R}  \cos(\theta) $'

    return p_phi_th, label_th




#--------------- Linera fit----------


def Lin_Func(x, a, b):
    return a*x + b


def Linear_Fit(x, y, params_hint):
    tab, pcov = curve_fit(Lin_Func, x, y, p0=params_hint)
    return tab, np.sqrt(np.diag(pcov))
    
    
















