#!/usr/bin/env python3
"""
Threshold behavior in the two-dimensional bistable equation

We leverage radial symmetry and illustrate threshold behavior in the
two-dimensional bistable equation given by

  1/xi (xi u')' + f(u) = 0

We use a planar system where the initial value is given by a small value of xi
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
def F(u, alpha):
	y = u*(1 - u)*(u - alpha)
	return y

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(xi, y, alpha):
	u, w = y
	du = -w/xi
	dw = xi*F(u, alpha)
	dy = np.array([du, dw])
	return dy

# =============================================================================
# DE Events
# =============================================================================
def de_event(t, y, alpha):
	u, w = y
	z = u*w
	return z
de_event.terminal = 1
de_event.direction = -1

# =============================================================================
# Main Simulation Function
# =============================================================================
def bistable_threshold_2D():
	a_list = [0.23, 0.2438355, 0.25]
	alpha = 0.1
	ximax = 100
	
	IVP_args = [alpha]
	fig1, ax1 = plt.subplots()
	ax1.set(xlabel = 'u', ylabel = 'w')
	fig2, ax2 = plt.subplots()
	ax2.set(xlabel=r'$\xi$', ylabel = r'u($\xi$)')
	
	for kc in range(len(a_list)):
		a = a_list[kc]
		xi0 = 0.01
		u0 = a - F(a, alpha)/4*xi0**2
		w0 = F(a, alpha)/2*xi0**2
		
		y0 = [u0, w0]
		soln = solve_ivp(de_rhs, [xi0, ximax], y0, args=IVP_args,
			events=de_event, dense_output=True)
		xif = max(soln.t)
		xi = np.linspace(xi0, xif, 2**8+1)
		U, W = soln.sol(xi)
		ax1.plot(U, W, label=f'a = {a:1.5f}')
		ax2.plot(xi, U, label=f'a = {a:1.5f}')
	ax1.legend(loc='upper right')
	ax2.legend(loc='upper right')
	plt.show()
	
		


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    bistable_threshold_2D()
