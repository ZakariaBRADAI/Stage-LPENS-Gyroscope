import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from IPython.display import HTML


def Gyro_Bloch(the, phi, psi):
    '''Sam's work. Unitary sphere, Sphère de Bloch'''

    x = np.sin(phi) * np.sin(the)
    y = -np.cos(phi) * np.sin(the)
    z = np.cos(the)

    return x, y, z


def Momentum_Weight(the, phi, psi, params):
    '''Moment associé au poids.'''
    g, m, h, _, _, _, _, _ = params
    momentum_the = m*g*h*np.sin(the)

    return momentum_the


def Momentum_Fe(t, the, phi, psi, params, forcing='X'):
    '''Moment associé à la force centrifuge.'''
    _, m, h, _, _, x0, Phi, omega = params
    match forcing:
        case 'X':
            momentum_the = m*h*x0*(omega**2)*np.cos(omega*t+Phi)* np.sin(phi)*np.cos(the)
            momentum_phi = m*h*x0*(omega**2)*np.cos(omega*t+Phi)* np.cos(phi)*np.sin(the)
        case 'XY':
            momentum_the = - m*h*x0*(omega**2)*np.sin(omega*t+Phi-phi)*np.cos(the)
            momentum_phi = m*h*x0*(omega**2)*np.cos(omega*t+Phi-phi)*np.sin(the)

    return momentum_the, momentum_phi


def Conjugated_Momentums(the, phi, psi, the_d, phi_d, psi_d, params):
    '''Moments conjugés.'''
    m, h, J1, J3 = params[1],  params[2], params[3],  params[4]
    J1_ = J1 + m*(h**2)
    p_psi = J3 * (psi_d + phi_d*np.cos(the))
    p_phi = J1_ * (np.sin(the)**2) * phi_d + np.cos(the)*p_psi
    p_the = J1_ * the_d
    return p_the, p_phi, p_psi


def p_psi_exp(params, CI):
    '''Calcule le P_psi initial.'''
    J3 = params[-4]
    p_psi0 = (CI[-1] + np.cos(CI[0]) * CI[-3]) * J3
    return p_psi0


def Get_Path(t, the, phi, psi, file_name="Gyro3D.png"):
    '''Génére la trajectoire du Gyro.'''

    x_t, y_t, z_t = Gyro_Bloch(the, phi, psi)
    
    i = 0
    f = len(t)

    layout = go.Layout(
        template="none",
        font=dict(family="Droid Serif"),
        scene=dict(
            xaxis_title=r"x",
            yaxis_title=r"y",
            zaxis_title=r"z",
            aspectratio=dict(x=1, y=1, z=1),
            camera_eye=dict(x=1.2, y=1.2, z=1.2),       
        ),
    )

    fig = go.Figure(layout=layout)
                    
    fig.add_scatter3d(x=[0], y=[0], z=[0], 
                      name="Contact point")
    
    fig.add_scatter3d(
        x=x_t[i:f],
        y=y_t[i:f],
        z=z_t[i:f],
        name="Path",
        mode="lines",
        line=dict(color=t, colorscale='Viridis', showscale=True),
    )

    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01))

    fig.update_layout(
        scene=dict(
            xaxis=dict(
                range=[-1, 1],
            ),
            yaxis=dict(
                range=[-1, 1],
            ),
            zaxis=dict(
                range=[-1, 1],
            ),
        )
    )
    fig.update_traces()
    fig.show()
    #fig.write_image(file_name)


def Compute_Carac_Pulsations(params, CI, forcing ='X'):
    '''Calcule les fréquences et pulsations théoriques attendues'''
    p_psi0 = p_psi_exp(params, CI)
    g, m, h, J1, J3, x0, p, omega_f = params
    #omega_f = 2 * np.pi * f

    omega_L_th = m * g * h / (p_psi0)
    match forcing:
        case 'X':
            omega_R_th = -0.5 * np.cos(p) * (m * h * x0 * (omega_f**2)) / p_psi0
        case 'XY':
            omega_R_th = -np.cos(p) * (m * h * x0 * (omega_f**2)) / p_psi0
    
    
    return omega_L_th, omega_R_th
    


def Gyro_Carac_Values(params, CI):
    '''Calcule les fréquences et pulsations théoriques attendues'''
    g, m, h, J1, J3, x0, p, f = params
    omega_L_th, omega_R_th = Compute_Carac_Pulsations(params, CI)
    print(f'Larmor Pulsation (th) : {omega_L_th : >+20_.3f}')
    print(f'Larmor Frequency (th) : {omega_L_th/(2*np.pi) : >+20_.3f}')
    print(f'Larmor Period (th) : {2 * np.pi / omega_L_th : >+20_.3f} \n')
        
    print(f'Rabi Pulsation (th) : {omega_R_th : >+20_.3f}')
    print(f'Rabi Period (th) : {2 * np.pi / omega_R_th : >+20_.3f}')
    print(f'Temps de montée (th) : {np.pi / omega_R_th : >+20_.3f} \n')
    
    
    rapport_freq = float(omega_L_th / omega_R_th)
    print(f'Rapport des pulsations Larmor/Rabi : {rapport_freq : >+20_.3f} \n')
        
    print(f'Rapport Approx Gyroscopique : {0.5 * J3 * ((2 * np.pi * CI[-1])**2) /  (m*g*h) : >+20_.3f} \n')
    
    return None


def Hamiltonian_Terms(t, the, phi, psi, the_d, phi_d, psi_d, params, forcing='Free'):
    '''Calcule les termes aparaissant dans l'Hamiltonien'''

    g, m, h, J1, J3, x0, p, w = params
    J1_ = J1 + m*(h**2)
    K = m * h * x0 * (w**2)
    p_the, p_phi, p_psi = Conjugated_Momentums(the, phi, psi, the_d, phi_d, psi_d, params)
    
    Epp = m*g*h*np.cos(the)
    Nutation = 0.5 * J1_ * (np.sin(the)**2) * (phi_d)**2
    Ec_theta = (p_the**2) / (2 * J1_)
    Ec_psi = (p_psi**2) / (2 * J3)
    match forcing:
        case 'Free':
            return Ec_theta, Ec_psi, Nutation, Epp
        case 'XY':
            E_ext = K*np.sin(the) * np.sin( w*t+p-phi )
            return Ec_theta, Ec_psi, Nutation, Epp, E_ext
        case 'X':
            E_ext = - K*np.sin(the) * np.sin(phi) * np.cos(w*t+p)
            return Ec_theta, Ec_psi, Nutation, Epp, E_ext


def Rabi_Ideal(t_burst, delta, Omega_R):
    '''Chevrons de Rabi idéaux.'''
    tab_Rabi_th = np.zeros((len(delta), len(t_burst)))
    for i in range(len(delta)):
        for j in range(len(t_burst)):
            pulsation = np.sqrt(delta[i]**2 + Omega_R**2) / 2
            num = (Omega_R**2) * np.sin( pulsation * t_burst[j] )**2
            den = Omega_R**2 + delta[i]**2
            tab_Rabi_th[i, j] = num / den
    return tab_Rabi_th


def Rabi_Freq_Modified(delta, Omega_R, Omega_L):
    return Omega_R * (1 + delta / Omega_L)**2


def Rabi_Assym(t_burst, delta, f_R, f_L):
    '''Chevrons de Rabi assymétriques.'''
    tab_Rabi_th = np.zeros((len(delta), len(t_burst)))
    for i in range(len(delta)):
        for j in range(len(t_burst)):
            Omega_R_mod = Rabi_Freq_Modified(delta[i], Omega_R, Omega_L)
            pulsation =  np.sqrt(delta[i]**2 + Omega_R_mod**2) / 2
            num = (Omega_R_mod**2) * np.cos( pulsation * t_burst[j] )**2
            den = Omega_R_mod**2 + delta[i]**2
            tab_Rabi_th[i, j] = num / den
    return tab_Rabi_th



'''
def Larmor_Freq_Modified_Moins(delta, f_L, t):
    return f_L - delta * np.sinc(delta * t / np.pi)


def Larmor_Freq_Modified_Plus(delta, f_L, t):
    return f_L+ delta * np.sinc(delta * t / np.pi)


def Rabi_Assym_DL(t_burst, delta, f_R, f_L):
    tab_Rabi_th = np.zeros((len(delta), len(t_burst)))
    for i in range(len(delta)):
        for j in range(len(t_burst)):
            f_L_mod = Larmor_Freq_Modified_Moins(delta[i], f_L, t_burst[j])
            f_R_mod = Rabi_Freq_Modified(delta[i], f_R, f_L_mod)
            pulsation = 2 * np.pi * np.sqrt(delta[i]**2 + f_R_mod**2) / 2
            num = (f_R_mod**2) * np.cos( pulsation * t_burst[j] )**2
            den = f_R_mod**2 + delta[i]**2
            tab_Rabi_th[i, j] = num / den
    return tab_Rabi_th
'''
