#!/usr/bin/env python3
"""
OU Process: We modify the Chapter 04 script single_particle_diffusion.py to
include drift due to deterministic velocity.

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
def xBrownian(D, kappa, dt, nt, Ntrials):
	# --- Simulation Tools
	X = np.zeros( (nt+1, Ntrials) )
	t = np.zeros(nt+1)

	for kt in range(nt):
		dX = -kappa*X[kt, :]*dt + np.sqrt(2*D*dt)*rng.normal(0, 1, Ntrials) 
		t[kt+1] = t[kt] + dt
		X[kt + 1, :] = X[kt, :] + dX
	return(t, X)

# =============================================================================
# Main Simulation Function
# =============================================================================
def OU_single_particle():
	# --- Parameters
	D = 1
	kappa = 1
	Np = 1000
	tf, Nt = 5, 2500
	
	# --- Get trajectories, do statistics
	t, X = xBrownian(D, kappa, tf/Nt, Nt, Np)
	Var = np.var(X, 1)
	
	# --- Do some plotting
	fig1, ax1 = plt.subplots()
	ax1.plot(t, X[:, :10])
	ax1.set(xlabel = 't', ylabel = 'Trajectories')
	
	fig2, ax2 = plt.subplots()
	ax2.plot(t, Var, label = 'Experiment')
	ax2.plot(t, D/kappa * (1 - np.exp(-4*kappa*D*t)), '--', label = 'Theory')
	ax2.set(xlabel = 't', ylabel = 'Variance')
	ax2.legend(loc = 'upper left')
	
	plt.show()
	


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    OU_single_particle()
