#!/usr/bin/env python3
"""
Phase Plane for the Steady State of the Fisher Equation

We find a numerical solution of the associated ODE, leveraging the event
switch in solve_ivp to stop integrating when needed. This simulation examines
three different cases: (1) Dirichlet BCs, (2) Robin BCs, and (3) BCs as a
result of a moving niche.

Produces: 
 - Figure 1: Phase portrait for ss Fisher Eqn as Figure 6.8
 - Figure 2: Phase portrait for ss Fisher Eqn as Figure 6.10
 - Figure 3: Phase portrait for ss Fisher Eqn as Figure 6.11
 - Figure 4: Survival/extinction curve as Figure 6.12

TO DO: Double check the simulation for moving niche BCs 
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# DE RHS (Dirichlet)
# =============================================================================
def de_rhs_D(t, y):
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
def D_BC_zero(t, y):
	v, w = y
	z = v
	return v

D_BC_zero.direction = -1	# Ensure v is decreasing 
D_BC_zero.terminal = True	# Trigger at first zero when z is decreasing

# =============================================================================
# DE RHS (Robin)
# =============================================================================
def de_rhs_R(t, y, delta):
	v, w = y
	dv = w
	dw = -v * (1 - v)
	dy = np.array([dv, dw])
	return dy

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
# DE RHS (Moving Niche)
# =============================================================================
def de_rhs_MN(t, y, c, lam_M):
	v, w = y
	dv = w
	dw =  - v * (1 - v)
	dy = np.array([dv, dw])
	return dy

# =============================================================================
# Moving Niche Boundary Conditions Event and Event Attributes
# =============================================================================
"""
The Moving Niche boundary conditions on the PDE give rise to a phase plane with
v > |w/lam_M| in traveling wave coordinates.  We ensure that the IVP solver returns 
such a solution.
"""
def MN_BC_zero(t, y, c, lam_M):
	v, w = y
	z = -lam_M * v + w
	return z

MN_BC_zero.direction = -1	# Ensure v is decreasing 
MN_BC_zero.terminal = True	# Trigger at first zero when z is decreasing


# =============================================================================
# Dirichlet BCs
# =============================================================================
def doDirichletBC(tMax):
	# --- Initialize Plot
	fig, ax = plt.subplots()
	ax.plot([0, 1], [0, 0], 'ok')
	ax.set(xlabel = 'v', ylabel = 'w')

	# --- Initialize IVP solve
	w0_max = np.sqrt(1/3)+2e-4
	w0_list = np.linspace(0.1, w0_max, 5)
	
	# --- Loop through initial conditions
	for jw0 in range(len(w0_list)):
		w0 = w0_list[jw0]
		Y_init = np.array([0, w0])
		soln = solve_ivp(de_rhs_D, [0, tMax], Y_init, 
						events = D_BC_zero, dense_output = True)
		tf = max(soln.t)
		t = np.linspace(0, tf, 2**8+1)
		y = soln.sol(t)
		v = y[0,:]
		w = y[1,:]
		ax.plot(v, w)
	ax.plot([0, 0], [-w0_max, w0_max], '--r')
	return fig, ax

# =============================================================================
# Robin BCs
# =============================================================================
def doRobinBC(tMax, delta):	

	# --- Initialize Plot
	fig, ax = plt.subplots()
	ax.plot([0, 1], [0, 0], 'ok')
	ax.set(xlabel = 'v', ylabel = 'w')

	# --- Initialize IVP solve
	w0_max = 0.5746745
	w0_list = np.linspace(0.1, w0_max, 5)
	ax.plot([0, w0_max/delta], [0, w0_max], '--r')
	ax.plot([0, w0_max/delta], [0, -w0_max], '--r')

	# --- Loop through initial conditions
	for jw0 in range(len(w0_list)):
		w0 = w0_list[jw0]
		Y_init = np.array([w0/delta, w0])
		soln = solve_ivp(de_rhs_R, [0, tMax], Y_init, args = [delta],
						events = R_BC_zero,  dense_output = True)
		tf = max(soln.t)
		t = np.linspace(0, tf, 2**8+1)
		y = soln.sol(t)
		v = y[0,:]
		w = y[1,:]
		ax.plot(v, w)
	return fig, ax

# =============================================================================
# Moving Niche
# =============================================================================
def doMovingNicheBC(tMax, delta, c):
	# --- Initialize Plot
	fig, ax = plt.subplots()
	ax.plot([0, 1], [0, 0], 'ok')
	ax.set(xlabel = 'v', ylabel = 'w')
	
	# --- Initialize IVP solve
	lam_P = (-c + np.sqrt(c**2 + 4*delta**2))/2
	lam_M = (-c - np.sqrt(c**2 + 4*delta**2))/2
	w0_max = 0.574375
	w0_list = np.linspace(0.1, w0_max, 5)
	ax.plot([0, w0_max/lam_P], [0, w0_max], '--r')
	ax.plot([0, -w0_max/lam_M], [0, -w0_max], '--r')

	# --- Loop through initial conditions
	for jw0 in range(len(w0_list)):
		w0 = w0_list[jw0]
		Y_init = np.array([w0/lam_P, w0])
		soln = solve_ivp(de_rhs_MN, [0, tMax], Y_init, args = [c, lam_M], 
						events = MN_BC_zero, dense_output = True)
		tf = max(soln.t)
		t = np.linspace(0, tf, 2**8+1)
		y = soln.sol(t)
		v = y[0,:]
		w = y[1,:]
		ax.plot(v, w)
	return fig, ax

# =============================================================================
# Critical value for survival
# =============================================================================
def doCritPlot(delta):
	# --- Variables
	c = np.linspace(0, 1.99, 2**8+1)
	y = 4/np.sqrt(4-c**2) * np.atan(np.sqrt( (4*delta**2 + c**2)/(4 - c**2)))
	y0 = 2*np.atan(delta)
	
	# --- Initialize plot
	fig, ax = plt.subplots()
	ax.set(xlim = (0, 2), ylim = (0, 4))
	ax.set(xlabel = 'C', ylabel = r'Y/Y$_0(\delta)$')
	
	# --- Annotate
	ax.annotate('survival', xytext = (0.5, 2), xy =(0.75, 2))
	ax.annotate('extinction', xytext = (1.25, 0.5), xy =(1.75, 0.5))
	ax.plot(c, y/y0)
	
	return fig, ax

# =============================================================================
# Main Simulation Function
# =============================================================================
def Fisher_ss_pp():
	# --- Parameters
	tmax = 20
	D = 0.1
	delta = 1/D
	c = 1
	
	# --- Dirichlet BCs
	fig1, ax1 = doDirichletBC(tmax)
	
	# --- Robin BCs
	fig2, ax2 = doRobinBC(tmax, delta)
	
	# --- Moving Niche
	fig3, ax3 = doMovingNicheBC(tmax, delta, c)
	
	# --- Critical Value
	fig4, ax4 = doCritPlot(np.sqrt(5))
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    Fisher_ss_pp()
