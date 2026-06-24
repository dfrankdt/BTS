#!/usr/bin/env python3
"""
Splitting Probability

We perform a random walk on the interval 0 < x < L, noting the distribution of 
leaving at x = L. We expect this distribution to be pi_L = x/L.

This script is based on splitting_probability.m

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()

# =============================================================================
# Stochastic Trajectory
# =============================================================================
def LR_exit(x0, L, D, dt):
	"""
	Step forward with
	
		x_t+1 = x_t + dx_t	

	where dx_t = sqrt(2*D*dt) N(0, 1) for some small dt.  Stop when x_t > L
	
	Inputs: 
		x0 (float): Initial position
		L (float): Width of boundary
		D (float): Diffusion coefficient
		dt (float): Temporal step 
		
	Uses:
		dx (float): Spatial step
		
	Outputs:
		x (float): 
	"""
	
	x = x0

	while (x < L)*(x > 0) > 0:
		dx = np.sqrt(2*D*dt)*rng.normal(0, 1)
		x = x + dx
	return x


# =============================================================================
# Main Simulation Function
# =============================================================================
def splitting_probability():
	# --- Simulation Parameters 
	D, L = 2, 1
	dt = 1/2**10
	Np = 1000
	xexit = np.zeros(Np)
	
	# --- Compute the distribution splitting probability
	nx0 = 2**4
	x0_stochastic = np.linspace(0.01, L - .01, nx0 + 1)
	pi_L = np.zeros(np.size(x0_stochastic))
	
	for kx0 in range(nx0+1):
		x0 = x0_stochastic[kx0]
		for kp in range(Np):
			xexit[kp] = LR_exit(x0, L, D, dt)
		pi_L[kx0] = np.sum(np.where(xexit>L, 1, 0))/Np
	
	# --- Do the plotting
	fig, ax = plt.subplots()
	ax.plot(x0_stochastic, pi_L, '.b', label = 'Stochastic Result')
	ax.plot([0, L], [0, 1], '--r', label='Deterministic Result')
	ax.set(xlabel = r'Initial Position, $x_0$', ylabel = r'Splitting Probability, $\pi_L$')
	ax.legend(loc='upper left')
	
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    splitting_probability()
