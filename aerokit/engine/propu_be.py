# coding: utf-8

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

# Calculs
def calculs(lambda_, pi_c, Tt4, pi_f, m):
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
    m_f = (lambda_ * m) / (lambda_ + 1)
    Tt17 = Tt2 * pi_f**((gamma - 1)/(gamma * eta_f))
    Pt17 = Pt2 * pi_f
    w_uf = m_f * cp * (Tt17 - Tt2)
    # Calculs CHP
    m_c = m / (lambda_ + 1)
    Tt3 = Tt2 * pi_c**((gamma - 1)/(gamma * eta_c))
    Pt3 = Pt2 * pi_c
    w_uc = m_c * cp * (Tt3 - Tt2)
    # Calculs chambre de combustion
    m_k = m_c * (((cp_ * Tt4) - (cp * Tt3))/((eta_comb * P_k) - (cp_ * Tt4)))
    alpha = m_k / m_c
    Pt4 = Pt3 * xi_cc
    # Calculs turbines
    Tt5 = Tt4 - ((w_uf + w_uc)/(m_c * (1 + alpha) * cp_ * eta_m))
    pi_t = (Tt4 / Tt5)**(gamma_ / (eta_turb * (gamma_ - 1)))
    Pt5 = Pt4 / pi_t
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
    Pt19 = Pt17 * xi_tuy
    Tt19 = Tt17
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
    print(w_uc)
    print(w_uf)
    print(Tt5)
    print("F = ",format(F,'.0f'))
    print("Fspe = ", format(Fspe,'.0f'))
    print("conso spe = ", format(mk_spe,'.10f'))
    print("eta_th = ",format(eta_th,'.3f'))
    print("eta_prop = ",format(eta_prop, '.3f'))
    print("eta = ",format(eta, '.3f'))
