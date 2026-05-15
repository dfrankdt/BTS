#!/usr/bin/env python3
"""
Travelling Waves of Fisher's Equation



"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(t, y, c):
	u, w = y
	du = w/c
	dw = -c * (w + u*(1 - u))
	dy = np.array([du, dw])
	return dy
	
# =============================================================================
# Solution Trajectory
# =============================================================================
def Usoln(tt, c, Uinit):
	tf = tt[-1]
	soln = solve_ivp(de_rhs, [0, tf], Uinit, args=[c], dense_output=True)
	Y = soln.sol(tt).T
	U = Y[:, 0]
	W = Y[:, 1]
	return U, W
	

# =============================================================================
# Main Simulation Function
# =============================================================================
def Fisher_waves_plot():
	# --- Global Parameters
	clist = np.array([0.5, 2, 3])
	Nt, tf = 2**8, 40
	
	# --- Simulation Vars
	t = np.linspace(0, tf, Nt+1)
	U0 = np.array([0.999, -0.001])	# Cannot start exactly at [1, 0] (Saddle)
	
	# --- Phase Plane
	fig1, ax1 = plt.subplots()
	x = np.linspace(-0.2, 1, 2**8+1)
	f = -x*(1-x)
	ax1.plot(x, f, '--k')
	ax1.set()
	ax1.set(xlabel = 'U', ylabel = 'W', title = 'U-W Phase Plane')

	# --- xi - U Plane
	fig2, ax2 = plt.subplots()
	ax2.set(xlabel = r'$\xi$', ylabel = r'$U_c(\xi)$', title='Traveling Wave Solution')
	for j in range(3):
		c = clist[j]
		# Phase plane for c
		u, w = Usoln(t, c, U0)
		ax1.plot(u, w)
		# Traveling wave for c
		traj_label = f'Trajectory with c = {c:1.2f}'
		ax2.plot(t - tf/2, u, label=traj_label)
	ax2.legend()
	plt.show()
	


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    Fisher_waves_plot()
