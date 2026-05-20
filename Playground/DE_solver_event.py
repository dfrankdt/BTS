#!/usr/bin/env python3
"""
DE Solver using events

This is a test script to see how the event switch works in ivp_solve.

I'm going to solve

  v'' + f(v) = 0,
  v(0) = 0, v'(0) = w0

the v - v' phase plane includes a homoclinic orbit through v = 1.  When the 
initial condition w0 lies inside of that homoclinic orbit, the solution is 
a periodic solution about the equilibrium solution v = 0.  I'd like to solve 
the corresponding differential equation while v > 0.

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
# Event and Event Attributes
# =============================================================================
def v_zero(t, y):
	v, w = y
	z = v
	return z

v_zero.terminal = True	# Trigger at first zero
v_zero.direction = 1	# Ensure v is decreasing (change to 1 to get the full cycle)

# =============================================================================
# Main Simulation Function
# =============================================================================
def DE_solve_event():
	fig, ax = plt.subplots()
	ax.plot([0, 1], [0, 0], 'ok')
	
	tmax = 20
	w0_list = np.linspace(0.1, 0.5, 8)
	for jw0 in range(len(w0_list)):
		w0 = w0_list[jw0]
		yinit = np.array([1e-10, w0])
		soln = solve_ivp(de_rhs, [0, tmax], yinit, 
				events = v_zero, dense_output = True)
		tf = max(soln.t)
		t = np.linspace(0, tf, 2**8+1)
		y = soln.sol(t)
		v = y[0,:]
		w = y[1,:]
		ax.plot(y[0,:], y[1,:])
		print(tf)
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    DE_solve_event()
