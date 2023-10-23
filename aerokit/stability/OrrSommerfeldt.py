
"""
    The ``OrrSommerfeld`` module
    =========================
 
    Provides 1D incompressible stability equations named Orr-Sommerfeld

    Orr-Sommerfeld equation by pseudo-spectral collocation method.

	\begin{equation*}
		\Bigg[ \frac{1}{i Re} (D^2 - \alpha^2)^2 - (\alpha U - \omega)(D^2 - \alpha^2) + \alpha U'' \Bigg] \hat{v} = 0 \qquad \text{Orr–Sommerfeld equation}
	\end{equation*}

	*   $D = \frac{d}{d y}$
	*   BCs: $\hat{v}(1) = \hat{v'}(1) = \hat{v}(-1) = \hat{v'}(-1) = 0$
	*   temporal theory: $\alpha \in \mathbb{R}$, $\omega \in \mathbb{C}$
	*   spatial theory: $\alpha \in \mathbb{C}$, $\omega \in \mathbb{R}$.
 
    :Example:
 
    >>> import aerokit.stability.OrrSommerfeld as OS
 
    Available functions
    -------------------
 
	.. note:: 
"""

import numpy as np
import aerokit.common.numspectral as ns

# def resol(alpha, omega, reynolds, DiffOp):
#     """
#     We write the Orr-Sommerfeld problem
#     in the following form:    mat.X = b

#     u   : base flow velocity
#     ddu : second derivative of u
#     e2  : second derivative matrix e2 = e*e
#     """

#     npts = DiffOp.npts
#     id_mat = np.eye(npts)
#     b   = np.zeros(npts, dtype=complex)
#     mat = np.zeros((npts, npts), dtype=complex)
#     ci = 1j
#     u = 1 - DiffOp.x**2
#     ddu = DiffOp.matder(2) @ u

#     for i in range(npts):
#       for j in range(npts):
#         mat[i,j] = alpha * (-u[i] * DiffOp.matder(2)[i,j] + (u[i] * alpha**2 + ddu[i]) * id_mat[i,j])

#     z1 = 1 / (ci * reynolds)
#     z2 = -2 * alpha**2 / (ci * reynolds) + omega
#     z3 = alpha**4 / (ci * reynolds) - omega * alpha**2
#     mat = mat + z1 * DiffOp.matder(4) + z2 * DiffOp.matder(2) + z3 * id_mat

#     mat[0,:] = DiffOp.matder(2)[0,:]
#     mat[1,:] = DiffOp.matder(1)[0,:]
#     mat[npts-2,:] = DiffOp.matder(1)[npts-1,:]
#     mat[npts-1,:] = 0
#     b[:] = 0
#     mat[-1,-1] = 1
#     b[0] = 1

#     sol_v = np.linalg.solve(mat, b)
#     csol = sol_v[0]
#     return csol, sol_v

# def spectre(alpha, reynolds, DiffOp, plot=True):
#     n = DiffOp.npts
#     mat1 = np.zeros((n, n), dtype=complex)
#     mat2 = np.zeros((n, n), dtype=complex)

#     ci = 1j
#     u = 1 - DiffOp.x**2
#     ddu = DiffOp.matder(2) @ u

#     id_matrix = np.eye(n) # identity matrix

#     for i in range(n):
#       for j in range(n):
#         mat1[i, j] = alpha * (-u[i] * DiffOp.matder(2)[i, j] + (u[i] * alpha**2 + ddu[i]) * id_matrix[i, j])

#     z1 = 1 / (ci * reynolds)
#     z2 = -2 * alpha**2 / (ci * reynolds)
#     z3 = alpha**4 / (ci * reynolds)
#     mat1 += z1 * DiffOp.matder(4) + z2 * DiffOp.matder(2) + z3 * id_matrix

#     z2 = -1
#     z3 = alpha**2
#     mat2 += z2 * DiffOp.matder(2) + z3 * id_matrix

#     # Boundary conditions
#     mat1[0, :] = 0
#     mat1[0, 0] = 1
#     mat2[0, :] = 0

#     mat1[1, :] = DiffOp.matder(1)[0, :]
#     mat2[1, :] = 0

#     mat1[n-2, :] = DiffOp.matder(1)[n-1, :]
#     mat2[n-2, :] = 0

#     mat1[n-1, :] = 0
#     mat1[n-1, n-1] = 1
#     mat2[n-1, :] = 0

#     l,v = eig(mat1, mat2)

#     if plot:
#         import matplotlib.pyplot as plt
#         plt.plot(l.real, l.imag, 'ob', markersize=6)
#         plt.hlines(y=0.0, xmin=-2, xmax=2, color='k', linestyle='--')
#         plt.xlabel('$\omega_r$')
#         plt.ylabel('$\omega_i$')
#         plt.gca().yaxis.label.set(rotation='horizontal', ha='right');
#         plt.title('Spectrum')
#         # plt.axis('square')
#         plt.axis([0.1, 1, -2, 0.5])
#         plt.show()
#     # else:
#     #     vp[:4] = l[:4]
#     return l, v



# ===============================================================
# automatic testing

if __name__ == "__main__":
    import doctest
    doctest.testmod()