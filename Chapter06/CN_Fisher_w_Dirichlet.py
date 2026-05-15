#!/usr/bin/env python3
"""
CN Scheme to solve Fisher's equation

We approximate the solution to Fisher's equation on 0 < x < L
(dimensionless) with zero Dirichlet boundary data using two
different values of L.

Notes
 - fix the initial profile
 - fine tune dt, dx 
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# =============================================================================
#  Nonlinearity
# =============================================================================
def F(u):
	y = u * (1-u)
	return y

# =============================================================================
#  CN Scheme
# =============================================================================
def doCN(u0, x, t):
	dx = x[1]-x[0]
	dt = t[1]-t[0]
	Nx = len(x) - 1
	Nt = len(t) - 1

	U = np.zeros( (Nx+1, Nt+1) )
	U[:,0] = u0

	# --- Matrices for performing CN
	gam = dt/(2*dx**2)
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)

	Acn = np.eye(Nx+1) - gam*D2
	Bcn = np.eye(Nx+1) + gam*D2
	
	Acn[0,:] = np.zeros(Nx+1)
	Acn[0, 0] = 1
	Acn[-1,:] = np.zeros(Nx+1)
	Acn[-1, -1] = 1

	# --- Initialization of CN method
	uk = u0
	for kt in range(Nt):
		y = Bcn@uk + dt*F(uk)
		y[0], y[-1] = 0, 0
		ukp1 = np.linalg.solve(Acn, y)
		uk = ukp1
		U[:,kt+1] = ukp1
	return U
# =============================================================================
#  Create Animation
# =============================================================================
def doMovie(x, t, U):
	Nt = len(t) - 1
	u0 = U[:,0]
        
	# Initialize movie
	fig, ax = plt.subplots()
	p_update = ax.plot([], [], 'b', label='Time Evolution')[0]
	p_init = ax.plot(x, u0, '--r', label='Initial Profile')
	ax.set(ylim=(0,1))
	ax.set(xlabel='x', ylabel='u(x, t)')
	ax.legend(loc='upper left')

	def update(frame):
		tk = t[frame]
		u = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(u)
		ax.set(title=f'Time t = {tk:.2f} s')
		ax.legend(loc='upper left')
		return(p_update)
        
	ani = manimation.FuncAnimation(fig=fig, func=update, 
	frames=range(0, Nt+1), interval=100)

	plt.show()

# =============================================================================
#  Main Simulation Function
# =============================================================================
def CN_Fisher_w_Dirichlet():
	# --- Global Parameter
	L = 5
	D = 1

	# -- Spatial and Temporal Scales
	Nt, Nx = 2**8, 2**8
	dt, dx = 0.002,  L/Nx
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, dt*Nt, Nt+1)

	# -- Initial profile and array for state variable
	u0 = np.tanh(x/(L/25)) * np.tanh(-(x - L)/(L/25))
	#u0 = (1 + np.tanh( (x - L/2)/(L/25) ))*(1 + np.tanh(-(x - L/2))/(L/25))/16
	#u0[0] = 0
	#u0[-1] = 0

	U = doCN(u0, x, t)

	doMovie(x, t, U)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
        CN_Fisher_w_Dirichlet()



