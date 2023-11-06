import math
import aerokit.aero.Isentropic as Is
import aerokit.aero.MassFlow as mf
import aerokit.common.defaultgas as defg

## Point fonctionnement Achile

Te = 290.   # Température extérieure
Pe = 101325. # Pression extérieure

T0 = 265. # Température totale jet

Mj = 1.4 # Nombre de Mach cible

# Calcul valeurs totales en entrée

defg.set_gamma(1.4) # useless since 1.4 is the default
r = 287.1

NPR = Is.PtPs_Mach(Mj)
P0 = NPR * Pe
rho0 = P0/(r*T0)

print('Q0:', T0, P0)

# Valeurs au col
Mc = 1.

# useless if you use (j) as reference but why not computing it
Tc = T0/Is.TtTs_Mach(Mc)
Pc = P0/Is.PtPs_Mach(Mc)
uc = Mc*math.sqrt(defg.gamma()*r*Tc)
rhoc = Pc/(r*Tc)

# Valeurs en entrée : on a besoin des surfaces en entrée et au col

Ri = 0.018 # m
Rc = 0.01122 # m

# isentropic relation for massflow conservation A/Sigma=cste
# we can directly compute it from Mj if you have Rj
Sigma_i = mf.Sigma_Mach(Mc)*(Ri/Rc)**2
Mi = mf.MachSub_Sigma(Sigma_i)


Ti = T0/Is.TtTs_Mach(Mi)
Pi = P0/Is.PtPs_Mach(Mi)

ui = Mi*math.sqrt(1.4*r*Ti)
rhoi = Pi/(r*Ti)
mi = math.pi*Ri**2*rhoi*ui

print('Qi:', ui, rhoi, Ti, Pi, mi)