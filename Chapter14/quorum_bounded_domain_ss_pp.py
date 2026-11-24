#!/usr/bin/env python3
"""
Quorum Sensing: We examine the first order system governing the steady-state
solution on the bounded domain.

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# Parameters
# =============================================================================
def doPars():
	# --- Parameters 
	s0 = 0.05
	gamma = 1.5
	alph = 0.0375
	K = 1
	d = 1
	D = 1
	return s0, K, gamma, alph, K, d, D
	
# =============================================================================
# Nonlinearity
# =============================================================================
def f(u):
	# --- Parameters
	s0, K, gamma, alph, K, d, D = doPars()
	y = s0 + gamma * u**2/(K**2 + u**2)
	return y

def fp(u):
	# --- Parameters
	s0, K, gamma, alph, K, d, D = doPars()
	yp = 2 * gamma * K * u/(K**2 + u**2)**2
	return yp

def fpp(u):
	# --- Parameters
	s0, K, gamma, alph, K, d, D = doPars()
	ypp = 2 * gamma * K**2 *(K**2 - 3*u**2)/(K**2 + u**2)**3
	return ypp

# =============================================================================
# u-w ODE RHS and Event
# =============================================================================
def uw_de_rhs(t, z, rho):
	# --- Parameters
	s0, K, gamma, alph, K, d, D = doPars()
	# --- State Variables
	u, w = z
	du = w/D
	dw = (fpp(u)*w**2/(d*D) + alph*u - (alph/d + rho/(1-rho))*f(u))
	dz = np.array([du, dw])

def uw_NC(t, z, rho):
	# --- Parameters
	s0, K, gamma, alph, K, d, D = doPars()
	# --- State Variables
	u, w = z
	w_NC = - np.sqrt(alph*D)*(u - f(u)/d)/(1 - fp(u)/d)
	return w - w_NC
uw_NC.terminal = 1
uw_NC.direction = -1	

# =============================================================================
# u-w Phase Plane
# =============================================================================
def get_uwPP(rho):
	# --- Parameters
	s0, K, gamma, alph, K, d, D = doPars()
	u = np.linspace(0, 3, 2**8 + 1)
	w = - np.sqrt(alph*D)*(u - f(u)/d)/(1 - fp(u)/d)
	
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.8, 6.4))
	ax1.plot(u, w, '--')
	ax1.set(xlabel = 'u', ylabel = r'u$_x$')
	
	u0_list = np.linspace(2, 3, 2**4+1)
	IVP_args = rho
	for ku0 in range(len(u0_list)):
		u0 = u0_list[ku0]
		z0 = [u0, 0]
		tmax = 100
		soln = solve_ivp(uw_de_rhs, [0, tmax], z0, args=[IVP_args], 
				events = uw_NC, dense_output = True)
		tf = max(soln.t)
		print(tf)
		t = np.linspace(0, tf, 2**8+1)
		u, w = soln.sol(t)
		ax1.plot(u, w)
	
	return fig, ax1, ax2

# =============================================================================
# Main Simulation Function
# =============================================================================
def quorum_bounded_domain_ss_pp():
	# --- Density
	rho = 0.04
	fig, ax1, ax2 = get_uwPP(rho)
	
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    quorum_bounded_domain_ss_pp()
