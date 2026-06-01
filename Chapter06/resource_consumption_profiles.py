#!/usr/bin/env python3
"""
Template

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(t, y, par_c, par_d):
	u, g, w = y
	c, d = par_c, par_d
	du = w
	dg = (c - w - c*(u + g))/d
	dw = -c*w -u*g 
	dy = np.array([du, dg, dw])
	return dy

# =============================================================================
# Main Simulation Function
# =============================================================================
def resource_consumption_profiles():
	# --- Global parameters
	c, d = 2, 20
	eps = 0.001
	lam = (-c+np.sqrt(c**2+4*d))/(2*d)
	
	# --- State variables and time span
	u0, g0, w0 = 1-eps, eps*lam*(c + lam), -eps*lam
	tf, Nt = 100, 100
	
	# --- Solution of ODE
	yinit = np.array([u0, g0, w0])
	soln = solve_ivp(de_rhs, [0, tf], yinit, args = [c, d], dense_output = True)
	
	# --- State variables solution
	t = np.linspace(0, tf, Nt+1)
	Y = soln.sol(t).T
	U = Y[:,0]
	G = Y[:,1]
	W = Y[:,2]
	
	# --- Plotting
	fig2, ax2 = plt.subplots()
	ax2.plot(t, U, 'b', label='U')
	ax2.plot(t, G, 'r', label='S')
	ax2.plot(t, W, 'y', label='W')
	ax2.plot(t, np.zeros(len(t)), '--c')
	ax2.set(xlabel = r'$\xi$') 
	ax2.legend(loc='upper right')

	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    resource_consumption_profiles()
