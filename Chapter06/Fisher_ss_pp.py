#!/usr/bin/env python3
"""
Phase Plane for the Steady State of the Fisher Equation

Note: Here we've leveraged the first integral to identify trajectories in the
case of Dirichlet boundary conditions, but this might change if the BCs are
Robin instead.  In the latter case we may need to approximate the solution of
the DE numerically.

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
def de_rhs(t, y):
	v, w = y
	dv = w
	dw = -v * (1 - v)
	dy = np.array([dv, dw])
	return dy

# =============================================================================
# Solution Trajectory
# =============================================================================
def Ysoln(tt, Uinit):
	tf = tt[-1]
	soln = solve_ivp(de_rhs, [0, tf], Uinit, dense_output=True)
	Y = soln.sol(tt).T
	V = Y[:, 0]
	W = Y[:, 1]
	return V, W
# =============================================================================
# First Integral
# =============================================================================
def F(v, w):
	z = 1/2*w**2 + 1/2*v**2 - 1/3*v**3
	return z


# =============================================================================
# Main Simulation Function
# =============================================================================
def Fisher_ss_pp():

	fig, ax = plt.subplots()
	ax.plot([0, 1], [0, 0], 'ok')
	ax.set(xlabel = 'v', ylabel = 'w')



	v = np.linspace(0, 1, 2**8+1)
	w = np.linspace(-0.6, 0.6, 2**8+1)
	V, W = np.meshgrid(v, w)
	
	w0_max = np.sqrt(1/3)
	w0_list = np.linspace(0.1, w0_max, 5)
	levels = 1/2*w0_list**2
	levels = np.sort(levels)
	
	ax.contour(V, W, F(V, W), levels)#, colors=['black'])
	
	plt.show()

	"""
	w0_list = np.linspace(0.1, 0.5, 5)
	t = np.linspace(0, 10, 2**8+1)
	fig, ax = plt.subplots()
	ax.plot([0, 1], [0, 0], 'ok')
	for jw0 in range(len(w0_list)):
		w0 = w0_list[jw0]
		Y_init = np.array([0, w0])
		v, w = Ysoln(t, Y_init)
		ax.plot(v, w, 'k')
	ax.set(xlim=(0, 1), ylim=(-0.5, 0.5))
	plt.show()
	"""

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    Fisher_ss_pp()
