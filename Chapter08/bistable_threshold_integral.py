#!/usr/bin/env python3
"""
Bistable Threshold Integral: 

We plot the relationship between alpha, the parameter of the nonlinearity in 
the bistable equation and Nc, an estimate of the stimulus required to produce
a traveling wave.

This integral is given by equation (8.23). We compute the value of the integral
by solving an initial value problem.

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
def de_rhs(t, y, c, alpha):
	u, w, v = y
	du = w
	dw = c*w - F(u, alpha)
	dv = u
	dy = np.array([du, dw, dv])
	return y

# =============================================================================
# DE Event
# =============================================================================
def de_event(t, y, c, alpha):
	u, w, v = y
	z = w*(1-u)
	return z
de_event.terminal = 1
de_event.direction = -1

# =============================================================================
# Main Simulation Function
# =============================================================================
def bistable_threshold_integral():
	# --- Parameters
	c = 0
	alpha_list = np.linspace(0.1, 0.5, 2**8)
	Int_value = np.zeros(alpha_list.shape)

	# --- ODE initialization
	u0 = 0.001
	tmax, Nt = 100, 2**8

	for kalpha in range(len(alpha_list)):
		alpha = alpha_list[kalpha]
		lam = (-c + np.sqrt(c**2 + 4*alpha))/2
		w0 = lam*u0
		y0 = [u0, w0, 0]
		IVP_args = [c, alpha]
		soln = solve_ivp(de_rhs, [0, tmax], y0, args=IVP_args,
			events=de_event, dense_output=True)

	U, V, IU = soln.y			
	fig, ax = plt.subplots()
	#ax.plot(alpha_list, Int_value)
	ax.plot(U, IU)
	ax.set(xlabel = r'$\alpha$', ylabel = r'$N_\alpha$')
	
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    bistable_threshold_integral()
    
    
