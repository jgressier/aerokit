
# coding: utf-8

# In[17]:


import math
import aerokit.aero.unsteady1D as uq
import aerokit.aero.riemann    as riem
import aerokit.aero.Isentropic as Is
import aerokit.aero.MassFlow   as mf
import aerokit.aero.degree     as deg
gam = 1.4


# Soit $\alpha$ l'angle max du distributeur, le rapport de section col est $A/A_c=1/\cos(\alpha)$. On calcule le nombre de Mach amont de cet écoulement bloqué

# In[13]:


alpha = 50.
AoAc  = 1./deg.cos(alpha)
M0    = mf.Mach_Sigma(AoAc, Mach=.2)
print("M0 = %4.3f pour A/Ac = %4.3f (%1.1f deg)"%(M0,AoAc,alpha))


# On prend pour référence une pression statique $p_0=1000\rm hPa$ au col et on définit un état initial 0 correspondant au Mach précédemment calculé.

# In[14]:


p0c = 1e5
Pt0 = p0c*Is.PtPs_Mach(Mach=1.)
#
p0 = Pt0/Is.PtPs_Mach(M0)
#
q0 = uq.unsteady_state(rho=1., u=M0*math.sqrt(gam*p0/1.), p=p0)
# check Mach and derived quantities
print("Q0: ",q0.Mach(), q0.Ptot(), q0.massflow())


# # hypothèse 1 : on impose l'état final
# 
# L'état final Q1f est défini par le ratio de Ptot et est supposé au même Mach (blocage) et même Ti

# In[15]:


Ptot_ratio=1.4
# Ptot variation with Ttot and Mach number constant : T & u csts
q1f = uq.unsteady_state(rho=q0.rho*Ptot_ratio, u=q0.u, p=q0.p*Ptot_ratio)
#
print("Q1f: ",q1f.Mach(), q1f.Ptot(), q1f.massflow())


# In[16]:


pb = riem.riemann_pb(q1f, q0)
print(pb.ustar(), pb.pstar())
print("post-Swave:", pb.qstarL().Ptot(), pb.qstarL().Mach(), pb.qstarL().massflow())
print("post-shock:", pb.qstarR().Ptot(), pb.qstarR().Mach(), pb.qstarR().massflow())


# # hypothèse 2 : on impose l'état à Ptot final mais Ps initiale
# 
# L'état final Q10 est défini par le ratio de Ptot mais garde la pression initiale

# In[21]:


# Ptot variation with Ps, Ttot constant : T & u csts
M10  = Is.Mach_PtPs(q0.Ptot()*Ptot_ratio/q0.p)
p10  = q0.Ptot()*Ptot_ratio/Is.PtPs_Mach(M10)
rho10= p10/(q0.p/q0.rho*Is.TtTs_Mach(M0)/Is.TtTs_Mach(M10))
q10 = uq.unsteady_state(rho=rho10, u=M10*math.sqrt(gam*p10/rho10), p=p10)
#
print("Q10: ",q10.Mach(), q10.Ptot(), q10.massflow())


# In[23]:


pb = riem.riemann_pb(q10, q0)
print(pb.ustar(), pb.pstar())
print("post-Swave:", pb.qstarL().Ptot(), pb.qstarL().Mach(), pb.qstarL().massflow())
print("post-shock:", pb.qstarR().Ptot(), pb.qstarR().Mach(), pb.qstarR().massflow())

