#!/usr/bin/env python3
"""
One Dimensional Run-and-Pause

This script is an algorithm to simulate movement of an object (say, bacterium)
that switches between moving in a one dimensional line with velocity v but 
randomly switches between left moving, stationary, and right moving. 

Figures produced:
	- Figure 1: Sample Trajectory
	- Figure 2: Mean Squared Displacement

Note that this code is based on one_d_run_and_pause.m
"""
import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()

# =============================================================================
# Run and Pause
# =============================================================================
def xtTrajectory(v, kon, koff, Tmax, Ntmax):
	"""
	Pass one particles through a direction switcher where a particle moves with 
	velocity v while moving and switches to moving with rate kon and to stationary
	with rate koff
	
	Inputs:
		v (float): Velocity particle moves while moving
		kon (float): Rate at which particle switches from stationary to moving
		koff (float): Rate at which particle switches from moving to stationary
		Tmax (float): Maximum time
		Ntmax (int): max number of time steps for the simulation
	Note that Ntmax is needed to populate the arrays t and x even though we do not 
	expect to use the entire array.
		
	Outputs:
		t (ndarray): time trajectory for a single particle
		x (ndarray): space trajectory for a single particle
	Note that time steps are random (via Gillespie) so t output is different each time
	"""
	x = np.zeros( Ntmax+1 )
	t = np.zeros( Ntmax+1 )
	
	kt = 0
	# -- Each cycle: start paused, wait time to switch, random to L/R, wait time to switch
	while t[kt] < Tmax:
		# -- One cycle
		R = rng.uniform(0, 1, 3)

		dt1 = -np.log(R[0])/kon				# Time to stay stationary
		pk = -1*(R[1]<1/2) + 1*(R[1]>1/2)	# Choose left or right
		dt2 = -np.log(R[2])/koff			# Time to move 

		dx1 = 0
		dx2 = pk*v*dt2

		x[kt+1] = x[kt] + dx1
		x[kt+2] = x[kt+1] + dx2
		t[kt+1] = t[kt] + dt1
		t[kt+2] = t[kt+1] + dt2
	
		kt = kt+2
	return t[:kt], x[:kt]		

# =============================================================================
# Main Simulation Function
# =============================================================================
def one_d_run_and_pause():
	"""
	Run through Np simulations, interpolate each onto uniform mesh for comparison
	"""
	
	# -- Parameters for the simulation -
	Np = 1000		# Particles
	kon = 7			# Rate to switch on
	koff = 8		# Rate to switch off
	v = 1			# Velocity when moving
	Tmax = 100		# Maximum time
	Ntmax = 5000	# Maximum steps (big, just to populate arrays)
	
	# -- Parameters for interpolating to uniform mesh -
	nt = 2**4
	tt = np.linspace(0, Tmax, nt+1)
	X = np.zeros( (nt+1, Np) )
	
	# -- Run through Np trajectories, compute the actual mean-squared displacement
	for kp in range(Np):
		t, x = xtTrajectory(v, kon, koff, Tmax, Ntmax)
		# -- Plot one such trajectory
		if kp == 0:
			plt.figure(1)
			plt.plot(t, x)
			plt.xlabel('time')
			plt.ylabel('x')
		# -- Do the interpolation
		X[:, kp] = np.interp(tt, t, x)

	# -- Do some statistics -
	Xms_actual = np.mean(X**2, 1)
	
	# -- We expect Deff = v^2/(koff) * (kon)/(kon + koff)
	Deff = v**2/(koff) * kon/(kon + koff)
	Xms_theory = 2*Deff*tt

	# -- Plot the Mean Squared Displacement
	plt.figure(2)
	plt.plot(tt, Xms_theory, '--r', label='Theoretical')
	plt.plot(tt, Xms_actual,'.k', label='Actual')
	plt.xlabel('time')
	plt.ylabel('Mean Squared Displacement')
	plt.legend()
	plt.show()
	
# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
	one_d_run_and_pause()

