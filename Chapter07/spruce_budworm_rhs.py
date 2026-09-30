#!/usr/bin/env python3
"""
Spruce Budworm

We plot the nonlinearities on right-hand side of the Spruce Budworm PDE (7.2)
to illustrate bistability of the system.

Produces:
 - Figure 1: Difference between growth and decay rates, fg - fd, as a function
   of the dimensionless parameter u/(K S), for varying sigma see Figure 7.1(a)
 - Figure 2: Division of kappa - sigma parameter space to illustrate bifurcation
   of number of zeros in fg - fd (Region I: one zero; Region II: three zeros),
   see Figure 7.1(b)
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Resultant
# =============================================================================
def R(x, y):
	R = 4*y**2*(y**2 + 1)**2*x**3
	R = R + 4*y**2*(3*y**2 - 5)*x**2
	R = R + (12*y**2 - 1)*x + 4
	return R

# =============================================================================
# Main Simulation Function
# =============================================================================
def spruce_budworm_rhs():
	# --- Parameters
	rB = 1.5
	K = 355
	beta = 43200
	alpha = 1.11
	rs = 0.1
	
	# --- Plotting: Bistability
	sigma_list = np.array([3, 4, 6, 8])
	a = alpha/K
	v = np.linspace(0, 1, 2**8+1)
	fig1, ax1 = plt.subplots()
	for sigma in sigma_list:
		fg = sigma*v*(1 - v)
		fd = v**2/(a**2 + v**2)
		f = fg - fd
		ax1.plot(v, f, label = rf'$\sigma = ${sigma:1.1f}')
	ax1.set(xlabel = r'$u/K_u S$', ylabel = r'$f_g(u) - f_d(u)$')
	ax1.plot([0, 1], [0, 0], '--k')
	ax1.legend()
		
	# --- Plotting: Resultant
	fig2, ax2 = plt.subplots()
	sigma = np.linspace(3, 6, 2**8+1)
	kappa = np.linspace(0, 0.2, 2**8+1)
	X, Y = np.meshgrid(sigma, kappa)
	ax2.contour(X, Y, R(X, Y), levels=[0])
	ax2.set(xlabel = r'$\sigma$', ylabel = r'$\kappa$')
	ax2.text(3.5, .07, 'I', size=18)
	ax2.text(4.5, .07, 'II', size=18)

	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    spruce_budworm_rhs()
