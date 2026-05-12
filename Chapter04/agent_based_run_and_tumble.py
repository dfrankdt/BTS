#!/usr/bin/env python3
"""
Agent based Run and Tumble

This script simulates run and tumble via an agent based model.  There
are Np particles that may be in one of three states, left moving, stationary
or right moving.

Figures produced:
	- Figure 1: Sample trajectories
	- Figure 2: Mean square displacement compared to theoretical

Note that this code is based on agent_based_run_and_tumble.m
"""

import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()

# =============================================================================
# Run and Tumble Cycle
# =============================================================================
def xRun_Tumble(kon, koff, v0, Np, Nt, dt):
	"""
	Pass Np particles through Nt timesteps
	
	Inputs:
		kon (float): Rate at which a particle switches to moving
		koff (float): Rate at which a particle switches to stationary
		v0 (float): velocity at which a particle moves when moving
		Np (int): Number of particles
		Nt (int): Number of timesteps
		dt (float): Length of timestep
		
	Uses:
		v (ndarray): Velocity vector
		cr (ndarray): Rates at which we leave a state
		s (int array): State vector (left, stationary, right)
		R (ndarray): Np uniformly distributed random numbers to determine switch
		pswitch (binary array): Identifies particles that switch
		ds (int array): Moves left and right to stationary, moves stationary left or right

	Outputs:
		t (ndarray): Times
		x (ndarray): Trajectories of Np particles
	"""
	v = v0*np.array([-1, 0, 1])
	cr = np.array([koff, kon, koff])

	t = np.zeros(Nt+1)
	x = np.zeros((Nt+1, Np))
	s = rng.integers(0, 2, Np, endpoint=True)

	for kt in range(Nt):
		# -- Move particles -
		dx = v[s]*dt
		t[kt+1] = t[kt] + dt
		x[kt+1,:] = x[kt, :] + dx
	
		# -- Determine which states switch -
		# -- Note: left moves +1, right moves -1, stationary moves left or right randomly
		R = rng.uniform(0, 1, Np)
		pswitch = (cr[s]*dt>R)
		ds = np.array([1, np.sign(rng.uniform(-1, 1)), -1], dtype=int)
		
		# -- Update state vector -
		s = s + ds[s]*pswitch
	return t, x

# =============================================================================
# Main Simulation Function
# =============================================================================
def agent_based_run_and_tumble():
	"""
	Identify the parameters needed, run the simulation, do some statistics and plotting
	"""

	# -- Parameters -
	Np = 2000
	Nt = 2000
	dt = 0.025
	v = .5
	kon = 2
	koff = 3

	t, x = xRun_Tumble(kon, koff, v, Np, Nt, dt)


	# -- Plot a few trajectories -
	plt.figure(1)
	plt.plot(t, x[:,:5])
	plt.xlabel('t')
	plt.ylabel('x')
	plt.title('Sample Trajectories')

	# -- Do some statistics -
	# -- We expect the effective diffusion coefficient Deff = v^2/(koff)*(kon/(kon + koff))
	Deff = v**2/(koff)*kon/(kon + koff)
	Xms_actual = np.mean(x**2, 1)
	Xms_theory = 2*Deff*t

	# -- Plot the theoretical mean-squared displacement and some of the calculated values
	plt.figure(2)
	plt.plot(t, Xms_theory, '--r', label='Theoretical')
	plt.plot(t[0:Nt+1:100], Xms_actual[0:Nt+1:100], '.k', label='Actual')
	plt.xlabel('t')
	plt.ylabel('Mean Square Displacement')
	plt.legend()

	plt.show()
	
# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    agent_based_run_and_tumble()

