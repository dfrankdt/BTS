#!/usr/bin/env python3
"""
Crank-Nicolson scheme to simulate the Bistable equation with a block, as described 
by equations (8.37) and (8.38). Below the simulation produces an animation.

We place a block of width Yb at location Xb. We find that that varying the width
of the block from 2.8 to 2.9 causes propagation failure.

Note: This script is based on CN_Bistable.m
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
def F(x, u, alpha, kappa, Xb, Yb):
	char_block = (x - Xb > 0)*(x - Xb < Yb)
	y_block = -kappa * u * char_block
	y_bistable = u * (1 - u) * (u - alpha) * (1 - char_block)
	y = y_block + y_bistable
	return y
# =============================================================================
# Crank-Nicolson Method
# =============================================================================
def doCN(x, t, uinit, alpha, kappa, X_block, Y_block):
	""" 
	Crank-Nicolson to simulate the Bistable Equation
		
		u_t = u_xx +  f(u)
		
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
		y = Bcn@uk + dt*F(x, uk, alpha, kappa, X_block, Y_block)
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
# Create Plots
# =============================================================================
def doPlot(x, t, U, ktskip):
	Nt = len(t) - 1
	uinit = U[:,0] 

	fig, ax = plt.subplots()
	ax.plot(x, uinit, '--r')
	ax.set(xlabel = r'$\xi$', ylabel = r'u($\xi$, $\tau$)')
	for kt in range(ktskip, Nt+1, ktskip):
		ax.plot(x, U[:, kt])
	plt.show()


# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_Bistable_block():
 
	# --- Parameters 
	L = 30			# Spatial Domain
	Tf = 200		# End time
	alpha = 0.25 	# Third zero of the cubic
	kappa = 0.1		# Decay constant for block
	Xb = 10			# Spatial location of block

	# --- Spatial and Temporal Scales 
	Nt, Nx = 2**8, 2**7
	dt, dx = Tf/Nt,  L/Nx
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, dt*Nt, Nt+1)

	# --- Initial profile for state variable 
	a = 0.42
	lam = 4.2
	u0_profile = a*(1/np.cosh((L - x)/lam))**2


	# --- Propagation: Solve, animate, graph
	Yb = 2.8
	U = doCN(x, t, u0_profile, alpha, kappa, Xb, Yb)
	doMovie(x, t, U, 2**3)
	doPlot(x, t, U, 2**4)

	# --- Propagation Failure: Solve, animate, graph
	Yb = 2.9
	U = doCN(x, t, u0_profile, alpha, kappa, Xb, Yb)
	doMovie(x, t, U, 2**3)
	doPlot(x, t, U, 2**4)

# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    CN_Bistable_block()





