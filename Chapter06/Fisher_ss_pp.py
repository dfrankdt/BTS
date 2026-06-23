#!/usr/bin/env python3
"""
Phase Plane for the Steady State of the Fisher Equation

We find a numerical solution of the associated ODE, leveraging the event
switch in solve_ivp to stop integrating when needed.

TO DO: Double check the trajectory in Figure 2
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
def de_rhs(t, y, delta):
	v, w = y
	dv = w
	dw = -v * (1 - v)
	dy = np.array([dv, dw])
	return dy

# =============================================================================
# Dirichlet Boundary Conditions Event and Event Attributes
# =============================================================================
"""
The Dirichlet boundary conditions on the PDE give rise to a phase plane with
v > 0 in traveling wave coordinates.  We ensure that the IVP solver returns 
such a solution.
"""
def D_BC_zero(t, y, delta):
	v, w = y
	z = v
	return z

D_BC_zero.direction = -1	# Ensure v is decreasing 
D_BC_zero.terminal = True	# Trigger at first zero when z is decreasing

# =============================================================================
# Robin Boundary Conditions Event and Event Attributes
# =============================================================================
"""
The Robin boundary conditions on the PDE give rise to a phase plane with
v > |w/delta| in traveling wave coordinates.  We ensure that the IVP solver returns 
such a solution.
"""
def R_BC_zero(t, y, delta):
	v, w = y
	z = delta * v + w
	return z

R_BC_zero.direction = -1	# Ensure v is decreasing 
R_BC_zero.terminal = True	# Trigger at first zero when z is decreasing

# =============================================================================
# Main Simulation Function
# =============================================================================
def Fisher_ss_pp():
	tmax = 20
	D = 0.1
	delta = 1/D 

	# --- Dirichlet Boundary Conditions
	fig1, ax1 = plt.subplots()
	ax1.plot([0, 1], [0, 0], 'ok')
	ax1.set(xlabel = 'v', ylabel = 'w')

	w0_max = np.sqrt(1/3)+2e-4
	w0_list = np.linspace(0.1, w0_max, 5)
	for jw0 in range(len(w0_list)):
		w0 = w0_list[jw0]
		Y_init = np.array([0, w0])
		soln = solve_ivp(de_rhs, [0, tmax], Y_init, args=[delta],
						events = D_BC_zero, dense_output = True)
		tf = max(soln.t)
		t = np.linspace(0, tf, 2**8+1)
		y = soln.sol(t)
		v = y[0,:]
		w = y[1,:]
		ax1.plot(v, w)
	ax1.plot([0, 0], [-w0_max, w0_max], '--r')
	
	# --- Robin Boundary Conditions
	fig2, ax2 = plt.subplots()
	ax2.plot([0, 1], [0, 0], 'ok')
	ax2.set(xlabel = 'v', ylabel = 'w')

	w0_max = 0.567125
	ax2.plot([0, w0_max/delta], [0, w0_max], '--r')
	ax2.plot([0, w0_max/delta], [0, -w0_max], '--r')
	w0_list = np.linspace(0.1, w0_max, 5)
	for jw0 in range(len(w0_list)):
		w0 = w0_list[jw0]
		Y_init = np.array([w0/delta, w0])
		soln = solve_ivp(de_rhs, [0, tmax], Y_init, args = [delta],
						events = R_BC_zero,  dense_output = True)
		tf = max(soln.t)
		t = np.linspace(0, tf, 2**8+1)
		y = soln.sol(t)
		v = y[0,:]
		w = y[1,:]
		ax2.plot(v, w)

	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    Fisher_ss_pp()
