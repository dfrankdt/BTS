#!/usr/bin/env python3
"""
DE Solver Test

This is a test script to see whether I even know how to use the DE solver

I'm going to solve a mass-spring equation of the form

  u'' + alpha^2 u = 0,
  u(0) = u0, u'(0) = 0

We know the solutions are of the form u(t) = u0 * cos(wt).  I'm interested
in using a DE event to stop the integration after one half a cycle, which 
should occur at tf = pi/w. 

The distance the mass travels in this time is u(pi/w) - u(0) = 2 u0.  We should
be able to determine this numerically by integrating v I guess.
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
def de_rhs(t, y, alpha):
	u, v, w = y
	du = v
	dv = -alpha**2 * u
	dw = u + v 
	dy = np.array([du, dv, dw])
	return dy

# =============================================================================
# Event and Event Attributes
# =============================================================================
def v_zero(t, y, alpha):
	u, v, w = y
	z = v
	return z

v_zero.terminal = 1 # Trigger at first zero
v_zero.direction = 1	# Ensure v is decreasing (change to 1 to get the full cycle)

# =============================================================================
# Main Simulation Function
# =============================================================================
def DE_solver():
	alpha = 3
	tmax = 20

	fig, ax = plt.subplots()
	
	u0_list = np.linspace(0.1, 0.5, 8)
	IVP_args = [alpha]
	for ju0 in range(len(u0_list)):
		u0 = u0_list[ju0]
		yinit = np.array([u0, 0, 0])
		soln = solve_ivp(de_rhs, [0, tmax], yinit, args = IVP_args,
				events = v_zero, dense_output = True)
		tf = max(soln.t)
		t = np.linspace(0, tf, 2**8+1)
		u, v, UI = soln.y
		print(soln.sol)
		ax.plot(u, v)
		print(UI)
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    DE_solver()
