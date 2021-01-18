import numpy as np
import hades.aero.Supersonic as sup
import matplotlib.pyplot as plt
import time

mach0 = np.linspace(1.2, 10., 201)
omeg  = sup.PrandtlMeyer_Mach(mach0)
start = time.clock()
mach1 = sup.old_Mach_PrandtlMeyer(omeg)
time1 = time.clock()-start
start = time.clock()
mach2 = sup.Mach_PrandtlMeyer(omeg)
time2 = time.clock()-start

print(time1, np.linalg.norm(mach1-mach0))
print(time2, np.linalg.norm(mach2-mach0))
plt.plot(mach0, mach1, mach0, mach2)
plt.show()