#!/usr/bin/env python3
"""
Diffusing particles with the possibility of degradation. Allow particles
to undergo Brownian motion. At each time step, allow for the possibility of
decay to a "dead" state, in which the particle no longer moves.

We collect the locations where particles "die".

Figures produced:
	- Figure 1: Histogram of location after decay with theoretical value
 
Note that this code is based on decay_probability.m
"""

import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()

# =============================================================================
# Brownian motion with decay
# =============================================================================
def diff_decay(D, alpha, dt, Np):
	"""
	When a particle moves, it does so according to a variance 2 D dt. At each time
	step, a particle may decay.
	
	Inputs:
		D (float): Diffusion coefficient
		alpha (float): Rate at which decay occurs
		dt (float): Time step
		Np (int): Number of particles
	
	Outputs: 
		X (ndarray): Array containing decay locations
	"""
	dx = np.sqrt(2*D*dt)
	X = np.zeros(Np)

	for j in range(Np):
		x = 0
		decay = 0
		while decay==0:
			x = x + dx*rng.normal(0, 1)
			if alpha*dt>rng.uniform(0, 1):
				decay = 1
		X[j] = x
	return X
# =============================================================================
# Main Simulation Function
# =============================================================================
def decay_probability():
	""" 
	Identify the parameters for the diffusion with decay simulation.  Run
	that simulation and do some statistics and plotting of the result.
	"""
	
	# -- Parameters -
	D = 1		# Diffusion coefficient
	alpha = 1	# Decay rate
	dt = 1e-4	# Time step
	Np = 1000	# Number of particles
	
	# -- Get distribution
	X = diff_decay(D, alpha, dt, Np)
	
	# -- Statistics -
	Xdist, bins = np.histogram(X, bins='auto', density=True)

	# -- Theoretical Distribution -
	x = np.linspace(-8, 8 , 2**8+1)
	p = np.sqrt(alpha/D)*np.exp(-np.sqrt(alpha/D)*np.abs(x))/2

	# -- Plotting -
	plt.figure(1)
	plt.stairs(Xdist, bins, fill=True)
	plt.plot(x, p, '--r')
	plt.xlabel('x')
	plt.ylabel('Distribution')
	plt.show()
	
# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    decay_probability()

