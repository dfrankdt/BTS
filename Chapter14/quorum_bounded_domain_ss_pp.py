#!/usr/bin/env python3
"""
Quorum Sensing: We examine the first order system governing the steady-state
solution on the bounded domain.

TO DO: Identify the outputs and better comment the code
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
	dw = (fpp(u)*w**2/(d*D) + alph*u - (alph/d + rho/(1-rho))*f(u))/(1 - fp(u)/d)
	dz = np.array([du, dw])
	return dz

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
	v = u - f(u)/d
	vx = (1 - fp(u)/d)*w
	
	fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(8.8, 6.4))
	ax1a.plot(u, w, '--')
	ax1a.set(xlabel = 'u', ylabel = r'u$_x$')
	ax1b.plot(v, vx, '--')
	ax1b.set(xlabel = 'v', ylabel = r'v$_x$')

	fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(8.8, 6.4))
	ax2a.set(xlabel = 'x', ylabel = 'u(x)')
	
	u0_list = np.linspace(0.1, 2.9, 2**8+1)
	uL_list = np.zeros(len(u0_list))
	L_list = np.zeros(len(u0_list))
	IVP_args = rho
	for ku0 in range(len(u0_list)):
		u0 = u0_list[ku0]
		z0 = [u0, 0]
		Lmax = 30
		soln = solve_ivp(uw_de_rhs, [0, Lmax], z0, args=[IVP_args], 
				events = uw_NC, dense_output = True)
		L = max(soln.t)
		x = np.linspace(0, L, 2**8+1)
		u, w = soln.sol(x)
		if np.mod(ku0, 2**3)==0:
			ax1a.plot(u, w)
			ax2a.plot(x, u)

			v = u - f(u)/d
			vx = (1 - fp(u)/d)*w
			ax1b.plot(v, vx)

		L_list[ku0] = L
		uL_list[ku0] = u[-1]

	ax2a.plot(L_list, uL_list, '--')
	ax2b.plot(L_list, u0_list, label = 'u(0)')
	ax2b.plot(L_list, uL_list, '--', label = 'u(L)')
	ax2b.legend(loc = 'upper right')	
	return fig1, ax1a, ax1b, fig2, ax2a, ax2b

# =============================================================================
# Main Simulation Function
# =============================================================================
def quorum_bounded_domain_ss_pp():
	# --- Density
	rho = 0.04
	fig1, ax1a, ax1b, fig2, ax2a, ax2b = get_uwPP(rho)
	
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    quorum_bounded_domain_ss_pp()
