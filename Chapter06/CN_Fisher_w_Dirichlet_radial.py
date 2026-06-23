#!/usr/bin/env python3
"""
CN Scheme to the 2D Fisher Equation (radial symmetry)

We solve the two-dimensional Fisher equation, leveraging radial symmetry, using
Crank-Nicolson

TO DO
 - Fix the second difference operator
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# Nonlinearity
# =============================================================================
def F(u):
	y = u*(1 - u)
	return y

# =============================================================================
# CN Solve
# =============================================================================
def doCN(r, t, uinit, D, UR):
	"""
	Crank-Nicolson to simulate the Fisher equation with radial symmetry
	
	  u_t = D  nabla^2 u + ru (K - U)
	  
	with u(R, t) = 0 and du/dr(0, t) = 0
	"""
	dr = r[1] - r[0]
	dt = t[1] - t[0]
	Nr = len(r) - 1
	Nt = len(t) - 1
	
	# --- Interior discretizations
	rp = (r[2:Nr+1] + r[1:Nr])/2
	rm = (r[1:Nr] + r[0:Nr-1])/2
	
	# --- Second derivative operator
	D2 = -2*np.eye(Nr+1)
	
	D2p = np.zeros(Nr)
	D2p[1:Nr] = rp/r[1:Nr]
	D2m = np.zeros(Nr)
	D2m[0:Nr-1] = rm/r[1:Nr]

	D2 = D2 + np.diag( D2p, k=1) + np.diag(D2m, k=-1)
	D2 = D2[0:Nr, 0:Nr]

	# --- Matrices for CN
	gam = D*dt/(dr**2)
	Acn = np.eye(Nr) - (gam/2)*D2
	Bcn = np.eye(Nr) + (gam/2)*D2
	
	# --- Boundary condition at r = R
	zbc = np.zeros(Nr)
	zbc[-1] = rp[-1]/r[-1]*UR
	zbc = gam*zbc

	# --- Initialization
	U = np.zeros( (Nr+1, Nt+1) )
	U[:, 0] = uinit
	
	# --- Steps
	uk = uinit[0:Nr]
	for kt in range(Nt):
		y = Bcn@uk + zbc + dt*F(uk)
		ukp1 = np.linalg.solve(Acn, y)
		U[0:Nr, kt+1] = ukp1
		U[-1, kt+1] = UR
		uk = ukp1
	return U
# =============================================================================
# Animation
# =============================================================================
def doMovie(r, t, U, ktskip):
	u0 = U[:,0]
	Nt = len(t) - 1

	# --- Initialization
	fig, ax = plt.subplots()
	p_init = ax.plot(r, u0, '--r', label = 'Initial Profile')
	p_update = ax.plot([], [], 'b', label = 'Time Evolution')[0]
	ax.set(xlabel = 'r', ylabel = 'u(r, t)', ylim=(0,1))
	ax.legend(loc = 'upper left')

	# --- Update
	def update(frame):
		tk = t[frame]
		u = U[:, frame]
		p_update.set_xdata(r)
		p_update.set_ydata(u)
		ax.set(title=f'Time t = {tk:.2f} s')
		return(p_update)
	ani = manimation.FuncAnimation(fig=fig, func=update, 
			frames=range(0, Nt+1, ktskip), interval=100)

	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_Fisher_w_Dirichlet_radial():
	# --- Global parameters
	R = 15 # or R = 5
	D = 1
	UR = 0
	
	# --- Discretization
	Nr = 2**6
	Nt = 2**8
	tf = 1
	r = np.linspace(0, R, Nr+1)
	t = np.linspace(0, tf, Nt+1)
	
	# --- Initial profile
	u0_profile = 0.2*np.exp(-r**2)
	U = doCN(r, t, u0_profile, D, UR)
	doMovie(r, t, U, 2**3)


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_Fisher_w_Dirichlet_radial()
