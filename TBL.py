import numpy as np
import matplotlib.pyplot as plt

def uplus_lin(yplus):
    kappa = .42
    E = 9.8
    A = (np.log(30*E)/kappa-5)/np.log(6)
    return yplus

def uplus_mix(yplus):
    kappa = .42
    E = 9.8
    A = (np.log(30*E)/kappa-5)/np.log(6)
    B = 5-A*np.log(5)
    return A*np.log(yplus)+B

def uplus_log(yplus):
    kappa = .42
    E = 9.8
    A = (np.log(30*E)/kappa-5)/np.log(6)
    return np.log(E*yplus)/kappa

fig, ax = plt.subplots(1, 1)

yp = np.linspace(1, 50., 100)
ax.plot(yp, uplus_lin(yp))
ax.plot(yp, uplus_mix(yp))
ax.plot(yp, uplus_log(yp))
#ax.set_xscale('log')
plt.show()

def uplus(yp):
    up = np.where(yp <5, uplus_lin(yp), uplus_mix(yp))
    up = np.where(yp <30, up, uplus_log(yp))
    return up

fig, ax = plt.subplots(1, 1)

yp = np.linspace(1, 100., 500)
ax.plot(yp, uplus(yp))
ax.set_xscale('log')
plt.show()

def utau_VF_estimate(u1, y1):
    kappa = .42
    E = 9.8
    return kappa*u1 / np.log(E*y1)
