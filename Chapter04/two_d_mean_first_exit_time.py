#!/usr/bin/env python3
"""
Two Dimensional First Exit Times

We extend the first_exit_times.py simulation to two dimensions, simulating 
Brownian motion through a stochastic differential equation and tracking the
first exit times from a certain domain, in this case a circle of radius R. 

Figures produced
 - Figure 1: A single realization of the stochastic differential equation
 - Figure 2: A histogram detailing the first exit times
 
Similar to first_exit_times.py, producing Figure 2 only requires returning the
exit time, so we do not compute the trajectory in its entirety. 


"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()

# =============================================================================
# Single Exit Time Trajectory
# =============================================================================
def two_d_exit_time_trajectory(R, D, dt):
	"""
	Step forward with
	
	  X_t+1 = X_t + dX_t
	  
	where dX_t = sqrt(2*D*dt) N(0, 1, (1, 2)) for some small dt. Stop when |X_t| > R
	
	Inputs:
		R (float): Radius of circle to define the exit
		D (float): Diffusion coefficient
		dt (float): Temporal step
	
	Uses:
		kt (int): Temporal counter
		dX (ndarray): Spatial step
		r (float): Distance from the origin
		Nt (int): Large number to prepopulate
		
	Returns:
		X (ndarray): Two dimensional trajectory
	"""
	r = 0
	kt, Nt = 0, int(10*R**2/(4*D)/dt)
	X = np.zeros( (2, Nt+1) )

	while r < R:
		dX = np.sqrt(2*D*dt)*rng.normal(0, 1, (1, 2))
		X[:, kt+1] = X[:, kt] + dX
		r = np.linalg.norm(X[:, kt+1])
		kt = kt + 1
	return X[0, :kt+1], X[1, :kt+1]

# =============================================================================
# Single Exit Time
# =============================================================================
def two_d_exit_time(R, D, dt):
	"""
	Step forward with
	
	  X_t+1 = X_t + dX_t
	  
	where dX_t = sqrt(2*D*dt) N(0, 1, (1, 2)) for some small dt. Stop when |X_t| > R
	
	Inputs:
		R (float): Radius of circle to define the exit
		D (float): Diffusion coefficient
		dt (float): Temporal step
	
	Uses:
		dX (ndarray): Spatial step
		r (float): Distance from the origin
		Nt (int): Large number to prepopulate
		
	Returns:
		t (float): Exit time
	"""
	r = 0
	t = 0
	X = np.zeros(2)
	while r < R:
		dX = np.sqrt(2*D*dt)*rng.normal(0, 1, (1, 2))
		X = X + dX
		t = t + dt
		r = np.linalg.norm(X)
	return t

# =============================================================================
# Main Simulation Function
# =============================================================================
def two_d_mean_first_exit_time():
	# --- Global Parameters
	Np = 2000
	R, D = 1, 1
	dt = 1/2**10

	# --- Plot one trajectory
	x, y = two_d_exit_time_trajectory(R, D, dt)
	t = np.linspace(0, 2*np.pi, 2**9+1)	
	fig1, ax1 = plt.subplots()
	ax1.plot(R*np.cos(t), R*np.sin(t), 'r')
	ax1.plot(x, y, 'b')
	ax1.set(xlabel = 'x', ylabel = 'y')
	ax1.set(title = r'Single Trajectory from $(x, y) = (0, 0)$')

	# --- Compute Np exit times without computing the trajectory
	texit = np.zeros(Np)
	for kp in range(Np):
		texit[kp] = two_d_exit_time(R, D, dt)

	# --- Do some statistics on the output (Figure 2)
	txdist, bins = np.histogram(texit, bins=40, density=False)
	fig2, ax2 = plt.subplots()
	ax2.stairs(txdist, bins, fill=True)
	ax2.set(xlabel = 't', title = r'First Exit Time from $(x, y) = (0, 0)$')

	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    two_d_mean_first_exit_time()
