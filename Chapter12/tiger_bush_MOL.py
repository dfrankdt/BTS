#!/usr/bin/env python3
"""
Method of Lines for Tiger Bush Model: We find the numerical solution to the
dimensionless Tiger Bush system (12.20) -- (12.21) using MOL.

Produces
 - Animation illustrating the traveling wave
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# PDE RHS
# =============================================================================
def pde_rhs(t, z, p):
	# --- Parameters and structure
	a, v, m, dx = p
	N = int( (len(z) - 2)/2 )
	w = z[:N+1]
	n = z[N+1:]
	
	# --- Nonlinearities and differences w/ periodic boundary conditions
	fw = a - w - w*n**2
	fn = w*n**2 - m*n
	
	wplus = np.roll(w, shift = -1)
	nplus = np.roll(n, shift = -1)
	nminus = np.roll(n, shift = 1)
	
	dwdt = fw + v*(wplus - w)/dx
	dndt = fn + (nplus - 2*n + nminus)/dx**2
	
	# --- Evaluate for dz/dt
	dz = np.zeros( (N+1) + (N+1) )
	dz[:N+1] = dwdt
	dz[N+1:] = dndt
	return dz

# =============================================================================
# Animation
# =============================================================================
def doMovie(x, t, U, W):
	# ----
	Nt = len(t) - 1
	uinit = U[:,0]
	winit = W[:,0]
	
	# --- Initialize
	fig, ax = plt.subplots()
	ax.set(xlim = (0, 8), ylim = (0, 10))
	ax.set(xlabel = r'$\xi$', ylabel = 'Density')
#	p_init_u = ax.plot(x, uinit, '--b')
#	p_init_w = ax.plot(x, winit, '--r')
	p_update_u = ax.plot([], [], 'b', label = 'Plant Density')[0]
	p_update_w = ax.plot([], [], 'r', label = 'Water Density')[0]
	ax.legend(loc = 'upper right')
	
	def update(frame):
		tk = t[frame]
		uk = U[:, frame]
		wk = W[:, frame]
		p_update_u.set_xdata(x)
		p_update_u.set_ydata(uk)
		p_update_w.set_xdata(x)
		p_update_w.set_ydata(wk)
		return(p_update_u, p_update_w)
	
	ani = manimation.FuncAnimation(fig = fig, func=update,
		frames=range(Nt+1), interval=100)
	
	return ani

# =============================================================================
# Main Simulation Function
# =============================================================================
def tiger_bush_MOL():
	# --- Parameters
	a = 7
	m = 2
	v = 20
	L = 7.85
	
	# --- Variables
	Nx = 2**8
	x = np.linspace(0, L, Nx+1)
	tf = 30
	
	# --- Initialization
	w0 = a*np.ones(Nx+1)
	n0 = m/w0 + 0.1*np.sin(2*np.pi*x/L)
	z0 = np.zeros( (Nx+1) + (Nx+1) )
	z0[:Nx+1] = w0
	z0[Nx+1:] = n0
	PDE_pars = [a, v, m, L/Nx]

	# --- MOL
	t = np.linspace(0, tf, 2**6+1)
	soln = solve_ivp(pde_rhs, [0, tf], z0, args=[PDE_pars], dense_output = True)
	Z = soln.sol(t)
	W = Z[:Nx+1, :]
	U = Z[Nx+1:, :]
	
	ani = doMovie(x, t, U, W)
	plt.show()	

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    tiger_bush_MOL()
