#!/usr/bin/env python3
"""
Cell Polarization: We examine ODE system describing the diffusion free model of cell polarization

Produces:
 - Figure 1: Right-hand side of the ODE, illustrating bistability, see Figure 12.10(a)
 - Figure 2: Steady states for different VT, see Figure 12.10(b)
 - Figure 3: Bifurcation diagram illustrating three regions of qualitatively
   different behavior, see Figure 12.11

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
def f(u, v, S, gamma):
	y = (S + gamma*u**2/(1 + u**2))*v - u
	return y

# =============================================================================
# RHS of DE - Illustrate bistability
# =============================================================================
def getRHS(gamma, VT, S_list):
	# --- Initiate indep variable, plot
	u = np.linspace(0, 1, 2**8+1)
	fig, ax = plt.subplots()

	# --- Step through S_list
	for kS in range(len(S_list)):
		S = S_list[kS]
		ax.plot(u, f(u, VT-u, S, gamma), label = f'S = {S:.2f}')
	ax.plot(u, np.zeros(len(u)), '--k')
	ax.set(xlabel = r'$u$', ylabel = r'$f(u, V_T - u)$', xlim = (0, 1), ylim = (-0.1, 0.1))
	ax.set(yticks = [-0.1, 0, 0.1])
	ax.legend(loc = 'upper right')
	return fig, ax

# =============================================================================
# Steady States
# =============================================================================
def getSS(gamma, VT_list):
	# --- Initiate indep variable, plot
	u = np.linspace(0, 1, 2**8+1)
	fig, ax = plt.subplots()

	# --- Step through VT_list
	for kV in range(len(VT_list)):
		VT = VT_list[kV]
		S = u/(VT - u) - gamma*u**2/(1 + u**2)
		ax.plot(S, u, label = rf'$V_T = {VT:1.1f}$')
	ax.set(xlabel = 'S', label = 'u', xlim = (0, 0.2))
	ax.set(xticks = [0, 0.1, 0.2])
	ax.legend(loc = 'upper right')
	
	return fig, ax

# =============================================================================
# Bifurcation Diagram
# =============================================================================
def getBifDiag():
	# --- Initiate indep variable, plot
	x = np.linspace(0.5, 4, 2**8+1)
	fig, ax = plt.subplots()
	
	# -- S(u) = S'(u) = 0 (note typo in text)
	y1 = 2/x*np.sqrt(x+1)
	
	# --- S'(u) = S''(u) = 0
	a = 27*x**2
	b = -2*(4*x**3 - 3*x**2 + 48*x + 32)
	c = 27*x**2
	
	z = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
	y2 = np.sqrt(z)
	
	# --- Plotting
	ax.plot(x, y1, x, y2)
	ax.set(xlabel = r'$\gamma$', ylabel = r'$V_T$', xlim=(0, 4), ylim=(0,4))
	
	return fig, ax
	
# =============================================================================
# Main Simulation Function
# =============================================================================
def cell_polarization_ode_system():
	# --- Parameters
	gamma = 3
	
	# --- RHS of DE
	VT = 1.3
	S_list = np.array([0.01, 0.05, 0.08])
	fig1, ax1 = getRHS(gamma, VT, S_list)

	# --- Steady states
	VT_list = np.array([1.1, 1.3, 1.4])
	fig2, ax2 = getSS(gamma, VT_list)

	# --- Bifurcation Diagram
	fig3, ax3 = getBifDiag()
	plt.show()
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	cell_polarization_ode_system()
