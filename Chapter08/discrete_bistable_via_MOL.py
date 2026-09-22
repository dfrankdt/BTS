#!/usr/bin/env python3
"""
Discrete Bistable via Method of Lines:

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# Nonlinearity
# =============================================================================
def f(x, a):
	# --- Cubic nonlinearity
	y = x*(1 - x)*(x - a)
	return y

# =============================================================================
# DE Structure
# =============================================================================
def de_rhs(t, u, pars):
	a, d = pars
	n = len(u)
	d2u = np.zeros(n)
	d2u[1:n-1] = d*(u[0:n-2] - 2*u[1:n-1] + u[2:n])
	d2u[0] = d*(-2*u[0] + u[1])
	d2u[-1] = d*(u[-2] - 2*u[-1])
	du = d2u + f(u, a)
	return du

# =============================================================================
# Plotting
# =============================================================================
def doPlot(t, u):
	fig, ax = plt.subplots()
	ax.plot(t, u)
	ax.set(xlabel = 't', ylabel = r'$u_n(t)$')
	return fig, ax	


# =============================================================================
# Main Simulation Function
# =============================================================================
def discrete_bistable_via_MOL():
	# --- Parameters
	alpha = 0.25
	n = 14
		
	# --- Solution structure: Propagation
	d = 0.02
	tend = 800
	u0 = np.zeros(n)
	u0[0:5] = 0.5
	IVP_args = [alpha, d]

	soln = solve_ivp(de_rhs, [0, tend], u0, args = [IVP_args], dense_output = True)
	t = np.linspace(0, tend, 2**8+1)
	u = soln.sol(t).T
	fig_wave, ax_wave = doPlot(t, u)

	# --- Solution structure: Failure
	d = 0.018
	tend = 300
	u0 = np.zeros(n)
	u0[0:5] = 0.5
	IVP_args = [alpha, d]

	soln = solve_ivp(de_rhs, [0, tend], u0, args = [IVP_args], dense_output = True)
	t = np.linspace(0, tend, 2**8+1)
	u = soln.sol(t).T
	fig_fail, ax_fail = doPlot(t, u)
	plt.show()
	
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    discrete_bistable_via_MOL()
