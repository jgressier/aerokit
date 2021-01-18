#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

import matplotlib.pyplot as plt
import numpy as np
import aerokit.aero.Isentropic as Is

# Grandeurs d'entrée
P0 = 227.00e2
T0 = 217.00
M0 = 0.8

# Paramètres
xi_e = 0.97
eta_c = 0.90
eta_f = 0.89
eta_comb = 0.99
xi_cc = 0.95
eta_m = 0.98
eta_turb = 0.88
xi_tuy = 0.97

P_k = 42800.00e3

# Gaz parfait
r = 287.00
r_ = 291.60
gamma = 1.40
gamma_ = 1.33

def calculs(lambda_, pi_c, Tt4, m, taux_meca):
    V0 = M0 * np.sqrt(gamma * r * T0)
    # Calculs grandeurs totales entrée
    Pt0 = P0 * (1 + ((gamma - 1) / 2) * M0**2)**(gamma/(gamma - 1))
    Tt0 = T0 * (1 + ((gamma - 1) / 2) * M0**2)
    # Calculs plan 2
    Tt2 = Tt0
    Pt2 = Pt0*xi_e
    # Gaz parfait
    cp = (r * gamma)/(gamma - 1)
    cp_ = (r_ * gamma_)/(gamma_ - 1)
    # Calculs CHP
    m_c = m / (lambda_ + 1)
    Tt3 = Tt2 * pi_c**((gamma - 1)/(gamma * eta_c))
    Pt3 = Pt2 * pi_c
    w_uc = m_c * cp * (Tt3 - Tt2)
    # Calculs chambre de combustion
    m_k = m_c * (((cp_ * Tt4) - (cp * Tt3))/((eta_comb * P_k) - (cp_ * Tt4)))
    alpha = m_k / m_c
    Pt4 = Pt3 * xi_cc
    # Calculs turbine compresseur
    Tt4_ = Tt4 - ((w_uc)/(m_c * (1 + alpha) * cp_ * eta_m))
    pi_tc = (Tt4 / Tt4_)**(gamma_ / (eta_turb * (gamma_ - 1)))
    Pt4_ = Pt4 / pi_tc
    # Calculs turbine fan
    P5 = P0
    T5 = Tt4_ * (P5 / Pt4_)**((gamma_ - 1) / gamma_)
    w_uf_theo = m_c * (1 + alpha) * cp_ * (T5 - Tt4_)
    Wu_f = w_uf_theo * taux_meca
    Tt5 = (Wu_f / (cp_ * m_c * (1 + alpha))) + Tt4_
    pi_tf = (Tt4_ / Tt5)**(gamma_ / (eta_turb * (gamma_ - 1)))
    Pt5 = Pt4_ / pi_tf
    # Calculs tuyère corps
    Tt9 = Tt5
    Pt9 = Pt5 * xi_tuy
    P9 = P0
    M9 = np.sqrt(((Pt9/P9)**((gamma_ - 1)/gamma_ ) - 1) * (2/(gamma_ - 1)))
    T9 = Tt9 / (1 + ((gamma_ - 1)/2)*M9**2)
    A9 = np.sqrt(gamma_ * r_ * T9)
    V9 = M9 * A9
    Fc = m_c * (V9 - V0)
    # Calculs fan
    m_f = (lambda_ * m) / (lambda_ + 1)
    Wf = -w_uf_theo * taux_meca * eta_m
    Tt13 = Wf/(m_f * cp) + Tt2
    pi_f = (Tt13 / Tt2)**(eta_f * gamma / (gamma - 1))
    Pt13 = Pt2 * pi_f
    # Calculs tuyère fan
    Pt19 = Pt13 * xi_tuy
    Tt19 = Tt13
    P19 = P0
    M19 = np.sqrt(((Pt19/P19)**((gamma - 1)/gamma ) - 1) * (2/(gamma - 1)))
    T19 = Tt19 / (1 + ((gamma - 1)/2)*M19**2)
    A19 = np.sqrt(gamma * r * T19)
    V19 = M19 * A19
    print("V9  inl", V9, M9)
    print("V19 inl", V19)
    Ff = m_f * (V19 - V0)
    # Calculs performances
    F = Fc + Ff
    Fspe = F/m
    mk_spe = (m_k / F)*10e6
    Qf = P_k * m_k
    delta_E = .5*(m_c*((1 + alpha)*V9**2 - V0**2) + m_f*(V19**2 - V0**2))
    eta_th = delta_E / Qf
    eta_prop = (F * V0) / delta_E
    eta = eta_th * eta_prop
    return Fspe, mk_spe, eta_th, eta_prop, eta, pi_f
    
def plot(X,Y,num_figure,title,x_label,y_label):
    figure = plt.figure(num=num_figure, figsize=(12,7))
    figure.suptitle(title)
    sub1 = plt.subplot(len(Y), 1, 1)
    plt.ylabel(y_label[0])
    plt.plot(X,Y[0])
    for i in range(1,len(Y)) :
        plt.subplot(len(Y),1,i+1, sharex=sub1)
        plt.plot(X,Y[i])
        plt.ylabel(y_label[i])
    plt.xlabel(x_label)
    return figure

def main(lambda_min,lambda_max,pi_c_min,pi_c_max,Tt4_min,Tt4_max):
    nb_valeurs = 20
    lambda_base = 6.1
    pi_c_base = 32.8
    Tt4_base = 1410
    m = 180
    taux_meca = 0.5
    y_label = ["Fspe","mk_spe","eta_th","eta_prop","eta_global"]
    lambda_array = np.linspace(lambda_min,lambda_max, num=nb_valeurs)
    pi_c_array = np.linspace(pi_c_min,pi_c_max, num=nb_valeurs)
    Tt4_array = np.linspace(Tt4_min,Tt4_max, num=nb_valeurs)
    taux_meca_array = np.linspace(0.2,0.8, num=nb_valeurs)
    results = [[0]*nb_valeurs for i in range(5)]
    for i in range(0,lambda_array.size) :
        result = calculs(lambda_array[i],pi_c_base,Tt4_base,m,taux_meca)
        results[0][i] = result[0]
        results[1][i] = result[1]
        results[2][i] = result[2]
        results[3][i] = result[3]
        results[4][i] = result[4]
    figure1 = plot(lambda_array,results,1,"Performances 1","Taux de dilution",y_label)
    for i in range(0,pi_c_array.size) :
        result = calculs(lambda_base,pi_c_array[i],Tt4_base,m,taux_meca)
        results[0][i] = result[0]
        results[1][i] = result[1]
        results[2][i] = result[2]
        results[3][i] = result[3]
        results[4][i] = result[4]
    figure2 = plot(pi_c_array,results,2,"Performances 2","Taux de compression coeur",y_label)
    for i in range(0,Tt4_array.size) :
        result = calculs(lambda_base,pi_c_base,Tt4_array[i],m,taux_meca)
        results[0][i] = result[0]
        results[1][i] = result[1]
        results[2][i] = result[2]
        results[3][i] = result[3]
        results[4][i] = result[4]
    figure3 = plot(Tt4_array,results,3,"Performances 3","Temperature sortie de chambre de combustion",y_label)
##    for i in range(0,taux_meca_array.size) :
##        result = calculs(lambda_base,pi_c_base,Tt4_base,m,taux_meca_array[i])
##        results[0][i] = result[0]
##        results[1][i] = result[1]
##        results[2][i] = result[2]
##        results[3][i] = result[3]
##        results[4][i] = result[4]
##    figure4 = plot(taux_meca_array,results,4,"Performances 4","Taux fan",y_label)
    figure1.show()
    figure2.show()
    figure3.show()
#    figure4.show()


def calculs_tau_optimum(lambda_, pi_c, Tt4, m):
    V0 = M0 * (gamma * r * T0)**(1/2)
    # Calculs grandeurs totales entrée
    Pt0 = P0 * (1 + ((gamma - 1) / 2) * M0**2)**(gamma/(gamma - 1))
    Tt0 = T0 * (1 + ((gamma - 1) / 2) * M0**2)
    # Calculs plan 2
    Tt2 = Tt0
    Pt2 = Pt0*xi_e
    # Gaz parfait
    cp = (r * gamma)/(gamma - 1)
    cp_ = (r_ * gamma_)/(gamma_ - 1)
    # Calculs fan
    tau_lambda = (cp_ * Tt4)/(cp * Tt0)
    tau_r = Tt0/T0
    tau_c = pi_c**((gamma - 1)/(eta_c * gamma))
    tau_f = (tau_lambda*tau_r*(tau_c - 1) - (tau_lambda/(tau_r*tau_c)) + lambda_*tau_r + 1)/(tau_r*(1 + lambda_))
    pi_f = tau_f**((eta_f*gamma)/(gamma - 1))
    Tt13 = Tt2 * tau_f
    Pt13 = Pt2 * pi_f
    m_f = (lambda_ * m) / (lambda_ + 1)
    wu_f = m_f * cp * (Tt13 - Tt2)
    # Calculs CHP
    m_c = m / (lambda_ + 1)
    Tt3 = Tt2 * pi_c**((gamma - 1)/(gamma * eta_c))
    Pt3 = Pt2 * pi_c
    wu_c = m_c * cp * (Tt3 - Tt2)
    # Calculs chambre de combustion
    m_k = m_c * (((cp_ * Tt4) - (cp * Tt3))/((eta_comb * P_k) - (cp_ * Tt4)))
    alpha = m_k / m_c
    Pt4 = Pt3 * xi_cc
    # Calculs turbine compresseur
    Tt5 = Tt4 - ((wu_c + wu_f)/(m_c * (1 + alpha) * cp_ * eta_m))
    pi_tc = (Tt4 / Tt5)**(gamma_ / (eta_turb * (gamma_ - 1)))
    Pt5 = Pt4 / pi_tc
    # Calculs tuyère corps
    Tt9 = Tt5
    Pt9 = Pt5 * xi_tuy
    P9 = P0
    M9 = (((Pt9/P9)**((gamma_ - 1)/gamma_ ) - 1) * (2/(gamma_ - 1)))**(1/2)
    T9 = Tt9 / (1 + ((gamma_ - 1)/2)*M9**2)
    A9 = (gamma_ * r_ * T9)**(1/2)
    V9 = M9 * A9
    Fc = m_c * (V9 - V0)
    # Calculs tuyère fan
    Pt19 = Pt13 * xi_tuy
    Tt19 = Tt13
    P19 = P0
    M19 = (((Pt19/P19)**((gamma - 1)/gamma ) - 1) * (2/(gamma - 1)))**(1/2)
    T19 = Tt19 / (1 + ((gamma - 1)/2)*M19**2)
    A19 = (gamma * r * T19)**(1/2)
    V19 = M19 * A19
    Ff = m_f * (V19 - V0)
    # Calculs performances
    F = Fc + Ff
    Fspe = F/m
    mk_spe = (m_k / F)*10e6
    Qf = P_k * m_k
    delta_E = (1/2)*(m_c*((1 + alpha)*V9**2 - V0**2) + m_f*(V19**2 - V0**2))
    eta_th = delta_E / Qf
    eta_prop = (F * V0) / delta_E
    eta = eta_th * eta_prop
    print(tau_f)
    print(tau_r)
    print(tau_c)
    print(tau_lambda)
    print(Tt13)
    #return [Fspe, mk_spe, eta_th, eta_prop, eta]



