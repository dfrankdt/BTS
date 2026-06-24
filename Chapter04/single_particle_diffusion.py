#!/usr/bin/env python3
"""
Single Particle Diffusion

We simulate the diffusion of a single particle via Brownian motion.

This script is based on single_particle_diffusion.m

"""

# =============================================================================
# Packages
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()

# =============================================================================
# Stochastic Differential Equation
# =============================================================================
def xBrownian(D, dt, nt, Ntrials):
	# --- Simulation Tools
	X = np.zeros( (nt+1, Ntrials) )
	t = np.zeros(nt+1)

	for kt in range(nt):
		dX = np.sqrt(2*D*dt)*rng.normal(0, 1, Ntrials) 
		t[kt+1] = t[kt] + dt
		X[kt + 1, :] = X[kt, :] + dX
	return(t, X)

# =============================================================================
# Main Simulation Function
# =============================================================================
def single_particle_diffusion():

	# --- Global Parameters
	D = 1/2
	dt = 1e-2
	Ntrials = 1000
	nt = 2500
	
	# --- Get the trajectories, do statistics
	t, X = xBrownian(D, dt, nt, Ntrials)
	var = np.var(X, 1)

	# --- Do some plotting
	fig1, ax1 = plt.subplots()
	ax1.plot(t, X[:,:9])
	ax1.set(xlabel = 'Time', ylabel = 'Sample Trajectories')

	fig2, ax2 = plt.subplots()
	ax2.plot(t, var, label='Experimental')
	ax2.plot(t, 2*D*t, '--', label='Theoretical')
	ax2.set(xlabel = 't', ylabel = 'Variance')
	ax2.legend(loc='upper left')

	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
	single_particle_diffusion()


