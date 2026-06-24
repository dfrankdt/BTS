#!/usr/bin/env python3
"""
Steady State Solutions

We plot steady state solutions of the diffusion equation with Robin conditions
represented by equation 5.24

This script is based on steady_state_solutions.m 
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def steady_state_solutions():
	# --- Global parameters
	Nx = 2**8
	x = np.linspace(0, 1, Nx+1)
	u0, uL = 1, 0
	Delta_list = np.array([0, 1, 10])

	# --- Do the plotting
	fig, ax = plt.subplots()
	for jDelta in range(np.size(Delta_list)):
		Delta = Delta_list[jDelta]
		u = 1/(1 + 2*Delta) * (uL - u0) * x + (Delta*(u0 + uL) + u0)/(1 + 2*Delta)
		ax.plot(x, u, label = rf'$\Delta = ${Delta:2d}')
	ax.plot([0, 1], [u0, uL], 'ok')
	ax.set(xlabel = r'$x/L$', ylabel = r'$u_s(x)$')
	ax.legend()
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    steady_state_solutions()
