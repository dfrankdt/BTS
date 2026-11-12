#!/usr/bin/env python3
"""
Seed Production: We construct solution profiles associated with the discretized
interpretation of the integrodifferential equation.

Produces
 - Figure 1: Snapshots of the solution u_n(x) 
 - Figure 2: Trace in the n-x plane where u_n(x) reaches a critical value
 - Figure 3: Animation of Figure 1
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# =============================================================================
# Nonlinearity
# =============================================================================
def f(u, A):
	y = u/(A + u)
	return y

# =============================================================================
# Production
# =============================================================================
def doProduction(x, u0, alpha, A, n):
	# --- Initialization
	dx = x[1] - x[0]
	Nx = len(x) - 1
	U = np.zeros( (Nx + 1, n + 1) )
	U[:, 0] = u0
	
	# --- Discretize Kernel
	y = x
	X, Y = np.meshgrid(x, y, indexing='ij')
	K = 1/2 * alpha * np.exp(- alpha * np.abs(X - Y) )
	
	# --- Iterate
	un = u0
	for iter in range(n):
		unp1 = K@f(un, A) * dx
		U[:, iter+1] = unp1
		un = unp1
	return U

# =============================================================================
# Snapshots
# =============================================================================
def doSnaps(x, U):
	fig, ax = plt.subplots()
	ax.set(xlabel = 'x', ylabel = r'$u_n(x)$', ylim = (0, 1.1*np.max(U)))
	ax.plot(x, U)
	
	return fig, ax

# =============================================================================
# Level Surfaces
# =============================================================================
def doLevels(x, U, ucrit):
	n = np.size(U, 1) - 1
	N, X = np.meshgrid(np.linspace(0, n, n+1), x)
	
	levels = [ucrit]
	fig, ax = plt.subplots()
	ax.contour(N, X, U, levels, colors=['black'], linestyles=['solid'])
	ax.set(xlabel = 'n', ylabel = 'x')
	return fig, ax

# =============================================================================
# Animation
# =============================================================================
def doMovie(x, U):
	# --- Initialization
	n = np.size(U, 1) - 1
	u0 = U[:, 0]
	
	fig, ax = plt.subplots()
	ax.set(xlabel = 'x', ylabel = r'$u_n(x)$', ylim = (0, 1.1*np.max(U)))
	p_init = ax.plot(x, u0, '--r')
	p_update = ax.plot([], [], '-b')[0]
	
	def update(frame):
		un = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(un)
		return(p_update)
		
	ani = manimation.FuncAnimation(fig = fig, func = update, 
		frames = range(n+1), interval = 100)
		
	return ani

# =============================================================================
# Main Simulation Function
# =============================================================================
def seed_production():
	# --- Parameters
	A = 0.5			# Hill Constant
	n = 35			# Number of iterations
	k, D = 1, 1 	# WLOG
	ulevel = 1 - A	# Find the the solution reaches this threshold
	alpha = np.sqrt(k/D)
	
	# --- Variables
	L, Nx = 100, 2**8
	dx = L/Nx
	x = np.linspace(0, L, Nx + 1)
	
	# --- Initialization and iteration
	u0 = 0.1 * 1 / np.cosh(0.5 * (x - L/2) )
	U = doProduction(x, u0, alpha, A, n)

	# --- Visualizations
	fig1, ax1 = doSnaps(x, U)
	fig2, ax2 = doLevels(x, U, ulevel)
	ani = doMovie(x, U)
	
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    seed_production()
