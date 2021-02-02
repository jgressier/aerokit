
# coding: utf-8

# <img align="right" src="./Downloads/cycle turbofan.tiff" width=400px>
# <h1 align="left" style="font-size:60px">Construction d'un cycle de Turbofan</h1>

# <h2>1. Importation des librairies python</h2>

# In[1]:


import numpy as np
import aerokit.aero.Isentropic as Is

# get_ipython().magic(u'matplotlib inline')
# get_ipython().magic(u'xmode verbose')

# ## 2. Introduction à la programmation objet python, définition de l'objet Gaz (gaz parfait)

# On propose de coder la classe $\textbf{Gaz}$, qui a pour arguments d'instance, $\gamma$, $r$, $P_k$ et une méthode get_cp(self).<br>
# Pour utiliser un objet en python il faut d'abords déclarer une classe comme en JAVA. 
# Au lieu de "$\textbf{def}$ function_name(args):" en python, on doit déclarer "$\textbf{class}$ class_name:"<br>
# Les attributs se déclarent comme des variable en python dans la class.<br>
# Pour déclarer une méthode on utilise "$\textbf{def}$ nom_de_methode(self, args):"<br>
# La méthode "$\textbf{__init__(self, args):}$" correspond au constructeur en JAVA.<br>
# L'argument "$\textbf{self}$" représente l'instance de l'objet (pointeur vers l'objet en memoire). On a donc 2 types d'arguments, les attributs qui sont communs à chaque objet (ex: M_0 = 0.8) et les attributs qui sont proprent à chaque objet (ex : self.gamma = 1.4).

# In[2]:


class Gaz:
    
    def __init__(self, gamma=1.4, r=287., P_k=0):
        self.gamma = gamma
        self.r = r
        self.P_k = P_k
        self.gsgmu = gamma / (gamma - 1.)
        self.gmusg = (gamma - 1.) / gamma
        
    def get_cp(self):
        return (self.r * self.gamma) / (self.gamma - 1.)


# ## 3. Construction du cycle de l'étage 1 à 4

# On construit maintenant une classe $\textbf{cycle_1_to_4}$ qui représentera un turbofan de l'étage 1 jusqu'à la sortie de chambre de combustion. Cette classe prend comme arguments : $P_0$, $T_0$, $M_0$, ainsi que les différents $\eta$ et $\xi$. Il y aura une méthode pour le calcule de chaque étage.<br><br>
# $\textbf{Rappels de cours :}$
# <ul>
#     <li>Étage 1 :</li>
#           <ul>
#                <li>Grandeurs totales : $$T_{i}=T\cdot \left(1 + \frac{\gamma -1}{2}\cdot M^2\right) \quad \& \quad P_{i}=P\cdot \left(1 + \frac{\gamma -1}{2}\cdot M^2\right)^{\frac{\gamma}{\gamma -1}}$$</li>
#           </ul>
#      <li>Étage 2 :</li>
#           <ul>
#                <li>Notion de perte de charge : $$P_{i2}=\xi \cdot P_{i1} \quad \& \quad T_{i2}=T_{i1}$$</li>
#           </ul>
#      <li>Étage 3 :</li>
#             <ul>
#                <li>Compression polytropique : $$\pi_c=\frac{P_{i3}}{P_{i2}}=\frac{T_{i3}}{T_{i2}}^{\frac{\gamma \cdot \eta_{pol}}{\gamma - 1}} \quad \& \quad \dot{w}_u=\dot{m}_c \cdot c_p \cdot \left(T_{i3} - T_{i2}\right)$$</li>
#           </ul>
#      <li>Étage 4:
#           <ul>
#                <li>Chambre de combustion : $$\alpha =\frac{\dot{m}_k}{\dot{m}} \quad \& \quad \dot{m}_c c^*_p\left(1+\alpha\right) T_{i4} - \dot{m}_c c_p T_{i3}=\eta_{comb} \dot{m}_k P_k \quad \& \quad P_{i4}=\xi \cdot P_{i3}$$</li>
#           </ul>
#      </li>
# </ul>

# In[3]:


class cycle_1_to_4:
    
    P0 = 227.00e2
    T0 = 217.00
    M0 = 0.8
    
    xi_e = 0.97
    eta_c = 0.90
    eta_f = 0.89
    eta_comb = 0.99
    xi_cc = 0.95
    eta_m = 0.98
    eta_turb = 0.88
    xi_tuy = 0.97
    
    def __init__(self, lambda_, pi_c, Tt4, m, g, g_fuel):
        self.lambda_ = lambda_
        self.pi_c = pi_c
        self.Tt4 = Tt4
        self.m = m
        self.g = g
        self.g_fuel = g_fuel
        self.current_stage_corps = 0
        
    def show_attributes(self):
        for attr in dir(self):
            if (not callable(getattr(self, attr)) and not attr.startswith('__')) :
                print("%s = %r" % (attr, getattr(self, attr)))
    
    def InitialValues(self):
        self.Pt0 = self.P0 * Is.PtPs_Mach(Mach=self.M0, gamma=self.g.gamma)
        self.Tt0 = self.T0 * Is.TtTs_Mach(Mach=self.M0, gamma=self.g.gamma)
        self.V0 = self.M0 * np.sqrt(self.g.gamma * self.g.r * self.T0)
        self.current_stage_corps = 1
        
    def stage_2(self):
        if self.current_stage_corps == 1:
            self.Tt2 = self.Tt0
            self.Pt2 = self.Pt0 * self.xi_e
            self.current_stage_corps = 2
        else:
            print ("Il faut initialiser les valeurs d'entrées à l'aide de la méthode InitialValues")
        
    def stage_3(self):
        if self.current_stage_corps == 2:
            self.m_c = self.m / (self.lambda_ + 1.)
            self.Tt3 = self.Tt2 * self.pi_c**(self.g.gmusg / self.eta_c)
            self.Pt3 = self.Pt2 * self.pi_c
            self.wu_c = self.m_c * self.g.get_cp() * (self.Tt3 - self.Tt2)
            self.current_stage_corps = 3
            
    def stage_4(self):
        if self.current_stage_corps == 3:
            self.m_k = self.m_c * (((self.g_fuel.get_cp() * self.Tt4) - (self.g.get_cp() * self.Tt3))/((self.eta_comb * self.g_fuel.P_k) - (self.g_fuel.get_cp() * self.Tt4)))
            self.alpha = self.m_k / self.m_c
            self.Pt4 = self.Pt3 * self.xi_cc
            self.current_stage_corps = 4
            
    def calculs_1_to_4(self):
        if self.current_stage_corps == 0:
            self.InitialValues()
            self.stage_2()
            self.stage_3()
            self.stage_4()
        if self.current_stage_corps == 1:
            self.stage_2()
            self.stage_3()
            self.stage_4()
        if self.current_stage_corps == 2:
            self.stage_3()
            self.stage_4()
        if self.current_stage_corps == 3:
            self.stage_4()
        


# ## 4. Construction du cycle l'étage 4 à la fin<br>
# On construit en suite une nouvelle classe $\textbf{cycle_5_to_9}$ qui hérite de la classe $\textbf{cycle_1_to_4}$ en déclarant "$\textbf{class}$ cycle_5_to_9(cycle_1_to_4):". Cela veut dire que la nouvelle classe héritera des attributs et des méthodes de la classe $\textbf{cycle_1_to_4}$.<br>
# Nous proposons 2 méthodes pour le calcul de la fin du cycle, une à $\pi_{fan}$ fixé et une autre à $\pi_{fan}$ calculé en renvoyant un pourcentage de la puissance récupérée dans la turbine BP vers le fan.<br><br>
# $\textbf{Rappels de cours :}$
# <ul>
#     <li>Turbines : $$\pi_t=\frac{P_{i4}}{P_{i5}}=\left(\frac{T_{i4}}{T_{i5}}\right)^{\frac{\gamma}{\eta_{pol}\left(\gamma - 1\right)}} \quad \& \quad \dot{w}_u=\dot{m}_c \left(1+\alpha \right) c^*_p \left(T_{i5} - T_{i4}\right)$$</li>
#     <li>Tuyères : $$T_{i9}\simeq T_{i5} \quad \& \quad P_{i9}=\xi \cdot P_{i5}$$</li>
#     <li>Poussée : $$F=\dot{m} \left(V_9 - V_0\right)$$</li>
#     <li>Puissance propulsive : $$P_F=F \cdot V_0$$</li>
#     <li>Variation d'énergie cinétique des gaz : $$\Delta E_c=\frac{1}{2} {m}_c \left(\left(1+\alpha\right)V^2_9 - V^2_0\right) + \frac{1}{2} {m}_f \left(V^2_{19} - V^2_{0}\right)$$</li>
#     <li>Puissance thermique de combustion : $$\dot{Q}_f=\dot{m}_k \cdot P_k$$</li>
#     <li>Rendements : $$\eta_T=\frac{\Delta E_c}{\dot{Q}_f} \quad \& \quad \eta_P=\frac{P_F}{\Delta E_c} \quad \& \quad \eta = \eta_T \cdot \eta_P$$</li>
#     <li>Grandeurs spécifique : $$F_{spe}=\frac{F}{\dot{m}}\quad \& \quad C_s=\frac{\dot{m}_k}{F}$$</li>
# </ul>

# ### 4.1. Méthode 1 : à $\pi_f$ fixé

# On déclare une classe qui hérite de la classe précedente

# In[4]:


class cycle_taux_fan_fixe(cycle_1_to_4):
    
    def __init__(self, lambda_, pi_c, Tt4, m, g, g_fuel, pi_f):
        self.lambda_ = lambda_
        self.pi_c = pi_c
        self.Tt4 = Tt4
        self.m = m
        self.g = g
        self.g_fuel = g_fuel
        self.pi_f = pi_f
        self.current_stage_corps = 0
        self.current_stage_fan = 0
        self.calculs_1_to_4()
        
    def stage_13(self):
        if self.current_stage_corps >= 2:
            self.m_f = self.m * self.lambda_ / (self.lambda_ + 1.)
            self.Tt13 = self.Tt2 * self.pi_f**(self.g.gmusg / self.eta_c)
            self.Pt13 = self.Pt2 * self.pi_f
            self.wu_f = self.m_f * self.g.get_cp() * (self.Tt13 - self.Tt2)
            self.current_stage_fan = 13
            
    def stage_19(self):
        if self.current_stage_fan == 13:
            self.Pt19 = self.Pt13 * self.xi_tuy
            self.Tt19 = self.Tt13
            self.P19 = self.P0
            self.M19 = Is.Mach_PtPs(self.Pt19/self.P19, self.g.gamma)
            self.V19 = Is.Velocity_MachTi(self.M19, self.Tt19, self.g.r, self.g.gamma)
            self.F_f = self.m_f * (self.V19 - self.V0)
            self.current_stage_fan = 19
            
    def stage_5(self):
        if self.current_stage_corps == 4 and self.current_stage_fan >= 13:
            self.Tt5 = self.Tt4 - ((self.wu_f + self.wu_c)/(self.m_c * (1. + self.alpha) * self.g_fuel.get_cp() * self.eta_m))
            self.pi_t = (self.Tt4 / self.Tt5)**(self.g_fuel.gsgmu / self.eta_turb)
            self.Pt5 = self.Pt4 / self.pi_t
            self.current_stage_corps = 5
            
    def stage_9(self):
        if self.current_stage_corps == 5:
            self.Tt9 = self.Tt5
            self.Pt9 = self.Pt5 * self.xi_tuy
            self.P9 = self.P0
            self.M9 = Is.Mach_PtPs(self.Pt9/self.P9, self.g_fuel.gamma)
            self.V9 = Is.Velocity_MachTi(self.M9, self.Tt9, self.g_fuel.r, self.g_fuel.gamma)
            self.F_c = self.m_c * (self.V9 - self.V0)
            self.current_stage_corps = 9
            
    def perfo(self):
        if (self.current_stage_corps == 9 and self.current_stage_fan == 19):
            self.F = self.F_f + self.F_c
            self.F_spe = self.F / self.m
            self.mk_spe = (self.m_k / self.F)*10e6
            Qf = self.g_fuel.P_k * self.m_k
            delta_E = .5*(self.m_c*((1. + self.alpha)*self.V9**2 - self.V0**2) + self.m_f*(self.V19**2 - self.V0**2))
            self.eta_th = delta_E / Qf
            self.eta_prop = (self.F * self.V0) / delta_E
            self.eta = self.eta_th * self.eta_prop
            
    def show_perfo(self):
        if self.current_stage_fan == 19 and self.current_stage_corps == 9:
            print("%s\t\t\t%.0f\t%s" % ("Poussée :", self.F, "N"))
            print("%s\t\t%.2f\t%s" % ("Poussée spécifique :", self.F_spe, "N.s/kg"))
            print("%s\t%.2f\t%s" % ("Consommation spécifique :", self.mk_spe, "g/kN/s"))
            print("%s\t\t%.3f\t%s" % ("Rendement thermique :", self.eta_th, ""))
            print("%s\t\t%.3f\t%s" % ("Rendement propulsif :", self.eta_prop, ""))
            print("%s\t\t%.3f\t%s" % ("Rendement global :", self.eta, ""))
    
    def calculs_5_to_9(self):
        if self.current_stage_fan == 0:
            self.stage_13()
            self.stage_19()
            self.stage_5()
            self.stage_9()
            self.perfo()
        if self.current_stage_fan == 13:
            self.stage_19()
            self.stage_5()
            self.stage_9()
            self.perfo()
        if self.current_stage_fan == 19:
            self.stage_5()
            self.stage_9()
            self.perfo()
        if self.current_stage_corps == 5:
            self.stage_9()
            self.perfo()
        if self.current_stage_fan == 19 and self.current_stage_corps == 9:
            self.perfo()
                               


# ### 4.1. Méthode 2 : à $\pi_f$ calculé

# In[9]:


class cycle_taux_fan_calcule(cycle_1_to_4):
    
    def __init__(self, lambda_, pi_c, Tt4, m, g, g_fuel, taux_meca=.5):
        self.lambda_ = lambda_
        self.pi_c = pi_c
        self.Tt4 = Tt4
        self.m = m
        self.g = g
        self.g_fuel = g_fuel
        self.taux_meca = taux_meca
        self.current_stage_corps = 0
        self.current_stage_fan = 0
        self.calculs_1_to_4()
        
    def stage_4_etoile(self):
        if self.current_stage_corps == 4:
            self.Tt4_etoile = self.Tt4 - ((self.wu_c)/(self.m_c * (1. + self.alpha) * self.g_fuel.get_cp() * self.eta_m))
            self.pi_t_c = (self.Tt4 / self.Tt4_etoile)**(self.g_fuel.gsgmu / self.eta_turb)
            self.Pt4_etoile = self.Pt4 / self.pi_t_c
            self.current_stage_corps = 5
            
    def stage_5(self):
        if self.current_stage_corps == 5:
            P5 = self.P0
            T5 = self.Tt4_etoile * (P5 / self.Pt4_etoile)**self.g_fuel.gmusg
            w_uf_theo = self.m_c * (1. + self.alpha) * self.g_fuel.get_cp() * (T5 - self.Tt4_etoile)
            self.wu_fan = w_uf_theo * self.taux_meca
            self.Tt5 = (self.wu_fan / (self.g_fuel.get_cp() * self.m_c * (1. + self.alpha))) + self.Tt4_etoile
            self.pi_t_f = (self.Tt4_etoile / self.Tt5)**(self.g_fuel.gsgmu / self.eta_turb)
            self.Pt5 = self.Pt4_etoile / self.pi_t_f
            self.current_stage_corps = 6
    
    def stage_9(self):
        if self.current_stage_corps == 6:
            self.Tt9 = self.Tt5
            self.Pt9 = self.Pt5 * self.xi_tuy
            self.P9 = self.P0
            self.M9 = Is.Mach_PtPs(self.Pt9/self.P9, self.g_fuel.gamma)
            self.V9 = Is.Velocity_MachTi(self.M9, self.Tt9, self.g_fuel.r, self.g_fuel.gamma)
            print("V9  cla", self.V9, self.M9)
            self.F_c = self.m_c * (self.V9 - self.V0)
            self.current_stage_corps = 9
            
    def stage_13(self):
        if self.current_stage_corps >= 6:
            self.m_f = (self.lambda_ * self.m) / (self.lambda_ + 1.)
            self.Tt13 = -self.wu_fan * self.eta_m /(self.m_f * self.g.get_cp()) + self.Tt2
            self.pi_f = (self.Tt13 / self.Tt2)**(self.eta_f * self.g.gsgmu)
            self.Pt13 = self.Pt2 * self.pi_f
            self.current_stage_fan = 13
            
    def stage_19(self):
        if self.current_stage_fan == 13:
            self.Pt19 = self.Pt13 * self.xi_tuy
            self.Tt19 = self.Tt13
            self.P19 = self.P0
            self.M19 = Is.Mach_PtPs(self.Pt19/self.P19, self.g.gamma)
            self.V19 = Is.Velocity_MachTi(self.M19, self.Tt19, self.g.r, self.g.gamma)
            print("V19 cla", self.V19)
            self.F_f = self.m_f * (self.V19 - self.V0)
            self.current_stage_fan = 19
            
    def perfo(self):
        if (self.current_stage_corps == 9 and self.current_stage_fan == 19):
            self.F = self.F_f + self.F_c
            self.F_spe = self.F / self.m
            self.mk_spe = (self.m_k / self.F)*10e6
            Qf = self.g_fuel.P_k * self.m_k
            delta_E = .5*(self.m_c*((1. + self.alpha)*self.V9**2 - self.V0**2) + self.m_f*(self.V19**2 - self.V0**2))
            self.eta_th = delta_E / Qf
            self.eta_prop = (self.F * self.V0) / delta_E
            self.eta = self.eta_th * self.eta_prop
            return self.F, self.F_spe, self.mk_spe, self.eta_th, self.eta_prop, self.eta
            
    def show_perfo(self):
        if self.current_stage_fan == 19 and self.current_stage_corps == 9:
            print("%s\t\t\t%.0f\t%s" % ("Poussée :", self.F, "N"))
            print("%s\t\t%.2f\t%s" % ("Poussée spécifique :", self.F_spe, "N.s/kg"))
            print("%s\t%.2f\t%s" % ("Consommation spécifique :", self.mk_spe, "g/kN/s"))
            print("%s\t\t%.3f\t%s" % ("Rendement thermique :", self.eta_th, ""))
            print("%s\t\t%.3f\t%s" % ("Rendement propulsif :", self.eta_prop, ""))
            print("%s\t\t%.3f\t%s" % ("Rendement global :", self.eta, ""))
            
    
    def calculs_5_to_9(self):
        if self.current_stage_corps == 4:
            self.stage_4_etoile()
            self.stage_5()
            self.stage_9()
            self.stage_13()
            self.stage_19()
            self.perfo()
        if self.current_stage_corps == 5:
            self.stage_5()
            self.stage_9()
            self.stage_13()
            self.stage_19()
            self.perfo()
        if self.current_stage_corps == 6:
            self.stage_9()
            self.stage_13()
            self.stage_19()
            self.perfo()
        if self.current_stage_corps == 9:
            self.stage_13()
            self.stage_19()
            self.perfo()
        if self.current_stage_fan == 13:
            self.stage_19()
            self.perfo()
        if self.current_stage_fan == 19 and self.current_stage_corps == 9:
            self.perfo()


# ## 5. Résultats

# In[13]:


# g = Gaz()
# g_fuel = Gaz(1.33,291.6,42800.e3)
# c1 = cycle_taux_fan_fixe(6.1,32.8,1410.,180.,g,g_fuel,1.6)
# c2 = cycle_taux_fan_calcule(11,41,1620.,230.,g,g_fuel,0.58)
# c1.calculs_5_to_9()
# c2.calculs_5_to_9()
# print("CFM56 :")
# c1.show_perfo()
# print "LEAP 1A :"
# c2.show_perfo()

