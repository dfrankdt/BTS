#!/usr/bin/env python3
"""
Crank-Nicolson scheme to simulate the Bistable equation. Below the simulation produces
an animation.

We note that the text uses initial data

 u0(x) = a sech^2(x/lam)

and cites a traveling wave with lam = 4.2 at a = 0.42 but no traveling wave at a = 0.41
We find that threshold to be slightly different, with no traveling wave at a = 0.40

Note: This script is based on CN_Fisher.m
"""

# =============================================================================
#  Packages
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation


# =============================================================================
#  Nonlinearity
# =============================================================================
def F(u, alpha):
	y = u * (1 - u) * (u - alpha)
	return y
# =============================================================================
# Crank-Nicolson Method
# =============================================================================
def doCN(x, t, uinit, alpha):
	""" 
	Crank-Nicolson to simulate the Bistable Equation
		
		u_t = u_xx +  u (1 - u) (u-a)
		
	with no-flux boundary conditions.
	"""
	dx = x[1]-x[0]
	dt = t[1]-t[0]
	Nx = np.size(x) - 1
	Nt = np.size(t) - 1

	U = np.zeros( (Nx+1, Nt+1) )
	U[:,0] = uinit

	# -- Matrices for performing CN
	gam = dt/(2*dx**2)
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	# -- No Flux BCs
	D2[0,0], D2[-1, -1] = -1, -1

	Acn = np.eye(Nx+1) - gam*D2
	Bcn = np.eye(Nx+1) + gam*D2
	
	# -- Initialization of CN method
	uk = uinit
	for kt in range(Nt):
		y = Bcn@uk + dt*F(uk, alpha)
		ukp1 = np.linalg.solve(Acn, y)
		U[:,kt+1] = ukp1
		uk = ukp1
	return U

# =============================================================================
# Create Movie
# =============================================================================
def doMovie(x, t, U, ktskip):
	# --- Initialize data structures
	Nt = np.size(t) - 1
	uinit = U[:,0]
	uMax = np.max(U)
	
	# --- Initialize movie
	fig, ax = plt.subplots()
	p_init = ax.plot(x, uinit, '--r', label='Initial Profile')
	p_update = ax.plot([], [], 'b', label='Time Evolution')[0]
	ax.set(ylim=(0, 1))
	ax.set(xlabel='x', ylabel='u(x, t)')
	ax.legend(loc='upper right')

	def update(frame):
	    tk = t[frame]
	    u = U[:, frame]
	    p_update.set_xdata(x)
	    p_update.set_ydata(u)
	    ax.set(title=f'Time t = {tk:.2f} s')
	    return(p_update)
        
	ani = manimation.FuncAnimation(fig=fig, func=update, 
			frames=range(0, Nt+1, ktskip), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_Bistable():
 
	# --- Parameters 
	L = 30			# Spatial Domain
	Tf = 40			# End time
	alpha = 0.25 	# Third zero of the cubic

	# --- Spatial and Temporal Scales 
	Nt, Nx = 2**10, 2**7
	dt, dx = Tf/Nt,  L/Nx
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, dt*Nt, Nt+1)

	# --- Initial profile for state variable (vary a to obtain a traveling wave)
	a = 0.41
#	a = 0.40
	lam = 4.2
	u0_profile = a*(1/np.cosh(x/lam))**2

	# --- Perform Crank-Nicolson
	U = doCN(x, t, u0_profile, alpha)

	# --- Create Movie 
	doMovie(x, t, U, 2**3)

# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    CN_Bistable()





