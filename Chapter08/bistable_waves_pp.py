#!/usr/bin/env python3
"""
Bistable Traveling Waves

We identify trajectories that yield traveling waves in the bistable equation
by examining the phase plane of the resulting ODE system.

TO DO
 - Work to obtain the saddle-saddle trajectory in the case c = c_crit
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
	u, w = y
	du = w
	dw = c*w - F(u, alpha)
	dy = np.array([du, dw])
	return dy

# =============================================================================
# DE Event
# =============================================================================
def de_event(t, y, c, alpha):
	# --- We keep our DE solution inside a box
	u, w = y
	z = w * (0.3 - w)
	return z
de_event.terminal = 1
de_event.direction = -1

# =============================================================================
# Main Simulation Function
# =============================================================================
def bistable_waves_pp():
	# --- Parameters
	alpha = 0.1
	c_list = [0, 0.56, 0.57, 1]
	
	# --- ODE initialization
	u0 = 0.001
	tmax, Nt = 100, 2**8
	
	# --- Plot initialization
	fig, ax = plt.subplots()
	ax.plot([0, 1], [0, 0], 'ok')
	ax.set(xlabel = 'U', ylabel = 'W')
#	ax.set(ylim=(-0.02, 0.2))
	
	
	# --- Simulations
	for kc in range(len(c_list)):
		c = c_list[kc]
		lam = (-c + np.sqrt(c**2 + 4*alpha))/2
		w0 = lam*u0
		y0 = [u0, w0]
		IVP_args = [c, alpha]
		soln = solve_ivp(de_rhs, [0, tmax], y0, args=IVP_args,
			events=de_event, dense_output=True)
		tz = max(soln.t)
		print(soln.t_events)
		t = np.linspace(0, tz, Nt+1)
		U, W = soln.sol(t)
		ax.plot(U, W, label = f'c = {c:1.2f}')
	
	# --- Critical value 
	c = np.sqrt(2)*(1/2 - alpha)
	c = 0.566662
	lam = (-c + np.sqrt(c**2 + 4*alpha))/2
	w0 = lam*u0
	y0 = [u0, w0]
	IVP_args = [c, alpha]
	soln = solve_ivp(de_rhs, [0, tmax], y0, args=IVP_args, 
		events = de_event, dense_output=True)
	tz = max(soln.t)
	t = np.linspace(0, tz, Nt+1)
	U, W = soln.sol(t)
	ax.plot(U, W, '--', label = f'c = c$^*$')

	ax.legend(loc='upper right')
	plt.show()
	
	# --- Trajectories in the xi plane
	fig, ax = plt.subplots()
	ax.plot(t, U, label=r'U($\xi$)')
	ax.plot(t, W, label=r'W($\xi$)')
	ax.set(xlabel=r'$\xi$')
	ax.legend()
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    bistable_waves_pp()
