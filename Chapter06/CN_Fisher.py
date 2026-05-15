#!/usr/bin/env python3
"""
Crank-Nicolson scheme to simulate the Fisher equation. Below the simulation produces
either a movie or a plot, depending on the comment.

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
def F(u, alpha, U0):
	y = alpha*u * (U0 - u)
	return y
# =============================================================================
# Crank-Nicolson Method
# =============================================================================
def doCN(uinit, x, t, Du, k, U0):
	""" 
	Crank-Nicolson to simulate the Fisher Equation
		
		u_t = Du u_xx + k u (U0 - u)
		
	with no-flux boundary conditions.
	"""
	dx = x[1]-x[0]
	dt = t[1]-t[0]
	Nx = np.size(x) - 1
	Nt = np.size(t) - 1

	U = np.zeros( (Nx+1, Nt+1) )
	U[:,0] = uinit

	# -- Matrices for performing CN
	gam = Du*dt/(2*dx**2)
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	# -- No Flux BCs
	D2[0,0], D2[-1, -1] = -1, -1

	Acn = np.eye(Nx+1) - gam*D2
	Bcn = np.eye(Nx+1) + gam*D2
	
	# -- Initialization of CN method
	uk = uinit
	for kt in range(Nt):
		y = Bcn@uk + dt*F(uk, k, U0)
		ukp1 = np.linalg.solve(Acn, y)
		uk = ukp1
		U[:,kt+1] = ukp1
	return U

# =============================================================================
# Create Movie
# =============================================================================
def doMovie(x, t, U):
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
	ax.legend(loc='upper left')

	def update(frame):
	    tk = t[frame]
	    u = U[:, frame]
	    p_update.set_xdata(x)
	    p_update.set_ydata(u)
	    ax.set(title=f'Time t = {tk:.2f} s')
	    return(p_update)
        
	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(0, Nt+1), interval=100)
	plt.show()

# =============================================================================
# Create Plots
# =============================================================================
def doPlots(x, t, U, Nskip):
	# --- Initialize data structures
	uinit = U[:,0]
	Nt = np.size(t) - 1
	
	# --- Refocus indices to choose just some of the t values
	t_indices = np.arange(0, Nt+1, Nskip)
	tk = t[ Nskip ]

	fig1, ax1 = plt.subplots()
	ax1.plot(x, U[:, t_indices])
	ax1.set(xlabel = 'x', ylabel = 'u(x, t)', 
			title = rf'Profiles every $\Delta t$ = {tk:1.2f} s')

	# --- Level set
	[T, X] = np.meshgrid(t, x)
	levels = [.5]
	fig2, ax2 = plt.subplots()
	ax2.contour(T, X, U, levels)

	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_Fisher():
 
	# --- Parameters 
	L = 100		# Spatial Domain
	Tf = 40		# End time
	Du = 1		# Diffusion coefficient
	alpha = 1	# Fisher rate constant 
	U0 = 1		# Carrying capacity

	# --- Spatial and Temporal Scales 
	Nt, Nx = 2**10, 2**7
	dt, dx = Tf/Nt,  L/Nx
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, dt*Nt, Nt+1)

	# --- Initial profile for state variable 
	uinit = 0.002*np.exp( -(x - L/2)**2 / (L/4) )

	# --- Perform Crank-Nicolson
	U = doCN(uinit, x, t, Du, alpha, U0)

	# --- Create Movie 
	# --- (comment out when producing a plot) 
	doMovie(x, t, U)
	
	# --- Create Plots 
	# --- (comment out when producing a movie) 
	# --- Nskip chooses profiles in time to plot from 0 to Nt by Nskip
	Nskip = 2**5
	doPlots(x, t, U, Nskip)

# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    CN_Fisher()





