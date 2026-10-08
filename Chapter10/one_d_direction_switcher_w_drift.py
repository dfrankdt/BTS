#!/usr/bin/env python3
"""
One Dimensional Direction Switcher with Drift: 
This script is an algorithm to simulate movement of an object (say, bacterium)
that switches between moving in a one dimensional line with velocity v but 
randomly switches direction by an exponential process with rate constants 
kplus and kminus.

Figures produced:
 - Figure 1: Sample Trajectories
 - Figure 2: Root Mean Squared Displacement (theoretical and actual)

We leverage the Chapter 04 one_d_direction_switcher.py to create this simulation

TO DO: The variance does not seem right
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()

# =============================================================================
# Direction Switcher
# =============================================================================
def xtTrajectory(v, kp, km, tmax, Ntmax):
	"""
	Pass one particle through a direction switcher where a particle moves with 
	velocity v while moving and switches with rates km and kp

	Inputs:
		v (float): Velocity particle moves while moving
		km (float): Rate at which particle switches from right to leftt motion
		kp (fl oat): Rate at which particle switches from left to right motion
		tmax (float): Maximum time
		Ntmax (int): Maximum Number of time steps for the simulation
                
	Outputs:
		t (ndarray): time trajectory for a single particle
		x (ndarray): space trajectory for a single particle
	Note that time steps are random (via Gillespie) so t output is different each time
	"""
	# --- Initialization
	x = np.zeros( Ntmax+1 )
	t = np.zeros( Ntmax+1 )
	kt = 0
	# --- Each cycle: move right, switch, move left, switch
	while t[kt] < tmax:
		# --- One cycle
		R = rng.uniform(0, 1, 2)

		dt1 = - np.log(R[0])/km # time moving right
		dt2 = - np.log(R[1])/kp # time moving left
		dx1 = v*dt1				# distance right	
		dx2 = -v*dt2			# distance left

		x[kt+1] = x[kt] + dx1
		x[kt+2] = x[kt+1] + dx2
		t[kt+1] = t[kt] + dt1
		t[kt+2] = t[kt+1] + dt2
	
		kt = kt+2
	return t[:kt], x[:kt]		


# =============================================================================
# Main Simulation Function
# =============================================================================
def one_d_direction_switcher_w_drift():
	# --- Global Parameters
	Np = 1000	    # particles
	kminus = 0.5    # switch rate R to L
	kplus = 1.5	    # switch rate L to R
	v = .5			# velocity while moving
	tmax = 100		# maximum time
	Ntmax = 5000	# maximum steps (big, just to preallocate arrays)
	
	# --- Parameters for interpolating to uniform mesh
	ntt = 2**4
	tt = np.linspace(0, tmax, ntt+1)
	X = np.zeros( (ntt+1, Np) )

	# --- Initialize Plot for ten trajectories
	fig1, ax1 = plt.subplots()
	ax1.set(xlabel = 'time', ylabel = 'x',
			title = 'One D Direction Switch with Drift')

	# -- Run through Np trajectories, compute the actual mean-squared displacement
	for kp in range(Np):
		t, x = xtTrajectory(v, kplus, kminus, tmax, Ntmax)
		# --- Plot ten such trajectories
		if kp < 10:
			ax1.plot(t, x)
		# --- Do the interpolation
		X[:, kp] = np.interp(tt, t, x)

	# --- Compute the mean, variance
	Xdrift = np.mean(X, 1)
	Xvar = np.var(X, 1)

	# --- Do some statistics
	Xrms = np.sqrt(np.mean(X**2, 1))

	
	# --- Drift: we expect veff = v*(kplus + kminus)/(kplus + kminus)
	veff = v*(kplus - kminus)/(kplus + kminus)
	Xdrift_theory = veff*tt
	fig2, ax2 = plt.subplots()
	ax2.plot(tt, Xdrift_theory, '--r', label='Theoretical')
	ax2.plot(tt, Xdrift,'.b', label='Actual')
	ax2.set(xlabel = 'time', ylabel = 'Drift',
				title = 'One D Direction Switch with Drift')
	ax2.legend()

	# --- Mean Squared Displacement: we expect Deff = v^2/(kplus + kminus), 
	Deff = v**2/(kplus + kminus)
	Xrms_theory = Deff*tt

	fig3, ax3 = plt.subplots()
	ax3.plot(tt, Xrms_theory , '--r', label='Theoretical')
	ax3.plot(tt, Xrms,'.b', label='Actual')
	ax3.set(xlabel = 'time', ylabel = 'Mean Squared Displacement',
				title = 'One D Direction Switch with Drift')
	ax3.legend()

	plt.show()
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    one_d_direction_switcher_w_drift()
