#!/usr/bin/env python3
"""
SIR Phase Portrait

Create the SIR model phase portrait shown in Figure 1.6

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# DE Right-Hand Side
# =============================================================================
def de_rhs(t, u, alpha, beta):
	"""
	Compute the right-hand side of the differential equation for the SIR model,
	  S' = -alpha SI
	  I' = alpha SI - beta I
	  R' = beta I
  	The constants alpha and beta are passed to the differential equation through
  	the parameter p. 
  	"""
	S, I, R = u
	du = np.array([-alpha*S*I, alpha*S*I - beta*I, beta*I])
	return du

# =============================================================================
# Subfigure (a)
# =============================================================================
def figure_a(tf, alpha, beta):
	Sinit_list = np.array([1.5, 2, 2.5])
	N_curves = np.size(Sinit_list)
	
	figa, axa = plt.subplots()
	for k in range(N_curves):
		Sinit = Sinit_list[k]
		uinit = np.array([Sinit, 0.001, 0])
		soln = solve_ivp(de_rhs, [0, tf], uinit, 
				args = [alpha, beta], dense_output=True)
		tt = np.linspace(0, tf, 2**9+1)
		uu = soln.sol(tt).T
		S = uu[:,0]
		I = uu[:,1]
		axa.plot(alpha*S/beta, alpha*I/beta)

	axa.plot([1, 1], [0, 1], '--k')
	# --- Add arrows
	axa.annotate("", xytext=(1.2, 0.8), xy = (0.8, 0.8),
				arrowprops = dict(arrowstyle="->"))

	axa.set(xlim = (0, 3), ylim = (0, 1))
	axa.set(xlabel = r'$\alpha s/\beta$', ylabel = r'$\alpha i/\beta$')

	
# =============================================================================
# Subfigure (b)
# =============================================================================
def figure_b():
	"""
	Plot (implicitly) the relationship between R0 = alpha s(0)/beta and s(infty)/s(0)
	"""
	
	# --- We cook up a variable x (representing s_infty/s0) that avoids 0 and 1
	Nx = 2**9
	dx = 1/Nx
	x = np.linspace(dx, 1-dx, Nx-1)
	R0 = np.log(x)/(x - 1)

	figb, axb = plt.subplots()
	axb.plot(R0, x)
	axb.set(xlabel = r'$R_0 = \alpha s(0)/\beta$', ylabel = r'$s(\infty)/s(0)$')

# =============================================================================
# Main Simulation Function
# =============================================================================
def SIR_pp():
	tf = 50
	alpha, beta = 1, 1
	figure_a(tf, alpha, beta)
	
	figure_b()
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    SIR_pp()
