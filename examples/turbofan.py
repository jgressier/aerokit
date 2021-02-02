import numpy                 as np
import aerokit.aero.Isentropic as Is
import matplotlib.pyplot     as plt
import aerokit.engine.turbofan as tf
#%matplotlib inline
plt.rcParams['figure.figsize'] = (10, 6)

M0=.5
T0=200.
P0=200e2
#
Ttmax = 1600.
Tt0   = T0*Is.TtTs_Mach(M0)
Tt3 = np.arange(1.2*Tt0, 0.5*Ttmax, 50.)
Opr = (Tt3/Tt0)**3.2
Bpr = np.arange(0.2, 10., .2)
AllOpr, AllBpr = np.meshgrid(Opr, Bpr)
#
model = tf.turbofan_adapt(AllOpr, Ttmax, AllBpr, AllBpr/(2.+AllBpr))
model.flight_conditions(T0, P0, M0)
model.update()
#
plt_etathp = plt.contour(AllOpr, model.spec_thrust(), model.thermoprop_efficiency(), 
						levels=np.arange(0, 1.,.05))
plt.clabel(plt_etathp, inline=True, fontsize=8)
#
plt_Bpr = plt.contour(AllOpr, model.spec_thrust(), AllBpr, linestyles='dashed', 
						levels=[0., 1., 2., 4., 8., 12.])
plt.clabel(plt_Bpr, inline=True, fontsize=8)
#
plt.show()