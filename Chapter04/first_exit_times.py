#!/usr/bin/env python3
"""
THIS DOES NOT LOOK RIGHT

First Exit Times 

This script simulates Brownian motion through a stochastic differential equation,
tracking the first exit times of a certain domain.  

Figures produced:
 - Figure 1: Trajectories of ten such motions
 - Figure 2: A histogram of first exit times
 - Figure 3: Comparison of mean first exit time to predicted quadratic
 
To produce Figure 1, we need the complete trajectory path. However the statistics 
and mean first exit time as a function of the nonzero initial position 0 < x0 < L
does not require the complete trajectory.
 
This script is based on first_exit_times.m

"""

import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()

# =============================================================================
# Single exit time trajectory
# =============================================================================
def exit_time_trajectory(x0, L, D, dt):
	"""
	Step forward with
	
		x_t+1 = x_t + dx_t	

	where dx_t = sqrt(2*D*dt) N(0, 1) for some small dt.  Stop when x_t > L
	
	Inputs: 
		x0 (float): Initial position
		dt (float): Time step
		L (float): Width of boundary
		D (float): Diffusion coefficient
	
	Uses:
		t0 (float): time constant for exit time
		k (integer): counter 
		t (float): curent time
		x (float): current position
		dx (float): space step
	
	Outputs:
		T (ndarray): Times
		X (ndarray): Trajectory
	"""
	t0 = L**2/(2*D)
	Nt = 10*int(t0/dt)	# Choose a large Nt so that dt*Nt >> t0

	t, x = 0, x0
	T, X = np.zeros(Nt), np.zeros(Nt)
	T[0], X[0] = 0, x0

	k = 0
	while x < L:
		# --- Step Forward
		dx = np.sqrt(2*D*dt)*rng.normal(0, 1)
		x = np.abs(x + dx)  # ensure reflecting boundary at x = 0
		t = t + dt
		
		# --- Collect Trajectories
		X[k+1] = x
		T[k+1] = t
		k = k + 1

	kexit = min(np.argwhere(X > L))	
	return T[0:k], X[0:k]

# =============================================================================
# Single exit time trajectory
# =============================================================================
def exit_time(x0, L, D, dt):
	"""
	Step forward with
	
		x_t+1 = x_t + dx_t	

	where dx_t = sqrt(2*D*dt) N(0, 1) for some small dt until x_t > L
	
	Inputs: 
		x0 (float): Initial position
		L (float): Width of boundary
		D (float): Diffusion coefficient
		dt (float): Time step
	
	Uses:
		x (float): Spatial position
		dx (float): Spatial step
		
	Outputs:
		t (float): Exit time
	"""
	t, x = 0, x0

	while x < L:
		# --- Step forward
		dx = np.sqrt(2*D*dt)*rng.normal(0, 1)
		x = np.abs(x + dx)  # ensure reflecting boundary at x = 0
		t += dt

	return t

# =============================================================================
# Main Simulation
# =============================================================================
def first_exit_times():

	# --- Simulation Parameters 
	D, L = 1, 1
	dt = 1/2**10
	Np = 1000
	texit = np.zeros(Np)

	# --- Plot ten trajectories (Figure 1)
	fig1, ax1 = plt.subplots()
	for kp in range(10):
		t, x = exit_time_trajectory(0, L, D, dt)
		ax1.plot(t, x)
		ax1.set(xlabel = 't', ylabel = 'x')
		ax1.set(title = 'Ten Trajectories')
		
	# --- Compute Np exit times without computing the trajectory 
	for kp in range(Np):
		texit[kp] = exit_time(0, L, D, dt)
	
	# --- Do some statistics on the output (Figure 2)
	txdist, bins = np.histogram(texit, bins=40, density=False)
	fig2, ax2 = plt.subplots()
	ax2.stairs(txdist, bins, fill=True)
	ax2.set(xlabel = 't', title = r'First Exit Time from $x=0$')
	
	# -- Collect exit time data for a number of initial points x0
	texit = np.zeros(Np)
	nx0 = 2**4
	x0_stochastic = np.linspace(0.01, L - .01, nx0 + 1)
	tm = np.zeros(nx0 + 1)
	for kx0 in range(nx0 + 1):
		for kp in range(Np):
			texit[kp] = exit_time(x0_stochastic[kx0], L, D, dt)
		tm[kx0] = np.mean(texit)
	
	# -- Plot the mean and compare it to the expected parabola (Figure 3)
	fig3, ax3 = plt.subplots()
	x0_deterministic = np.linspace(0, L, nx0 + 1)
	tm_deterministic = L**2/(2*D)*(1 - x0_deterministic**2/L**2)
	ax3.plot(x0_deterministic, tm_deterministic, '--r', label='Deterministic Result')
	ax3.plot(x0_stochastic, tm, '.b', label='Stochastic Result')
	ax3.set(ylim=(0,L**2/D))
	ax3.set(xlabel = 'Initial Position', ylabel = 'Mean Exit Time')
	ax3.legend(loc='upper left')

	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    first_exit_times()



