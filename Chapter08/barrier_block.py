#!/usr/bin/env python3
"""
Barrier Block:

Figures produced:
 - Figure 1: Trajectory (U-W plane), see Figure 8.6(a)
 - Figure 2: Trajectory (xi-U plane), see Figure 8.6(b)
 - Figure 3: Length of curve Y as a function of U(0), see Figure 8.7(a)
 - Figure 4: Critical blocking

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# Nonlinearities
# =============================================================================
def f(x, a):
	# --- Typical cubic nonlinearity
	y = x*(1 - x)*(x - a)
	return y

def F(x, a):
	# --- Integral of the cubic nonlinearity, zero of F
	y = -(x**4/4 - (a+1)/3*x**3 + a/2*x**2)
	x0 = 2*(a+1)/3 - np.sqrt(2)/3*np.sqrt( (1-2*a)*(2-a) )
	return x0, y

# =============================================================================
# DE Structure
# =============================================================================
def de_rhs_NL(x, z, a):
	u, w = z
	du = w
	dw = -f(u, a)
	dz = np.array([du, dw])
	return dz

def w_zero(x, z, a):
	u, w = z
	return w
w_zero.terminal = True
w_zero.direction = -1

# =============================================================================
# Figure 8.6 (a)
# =============================================================================
def getFig86a(U0):
	z0 = np.array([0.01, 0.01])
	soln = solve_ivp(de_rhs_NL, [-10, 10], z0, events = w_zero, dense_output = True)
	
	fig, ax = plt.subplots()
	ax.plot(soln.

# =============================================================================
# Main Simulation Function
# =============================================================================
def barrier_block():
	# --- Parameters
	


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    barrier_block()
