#!/usr/bin/env python3
"""
Agent based 2-dimensional SIR

This script simulates agent based process in which particles move randomly
but may interact according to S -> I with rate alpha

Figures produced:
	- Figure 1: Sample trajectories

This script is based on agent_SIR_1d.py
"""
# =============================================================================
# Packages
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation
rng = np.random.default_rng()

# =============================================================================
# React and Diffuse Cycle
# =============================================================================
def xyReactDiff(alpha, D, L, S0, Nt, dt):
	"""
	Pass Np = (L * S0) (must be even) particles through Nt timesteps. Particles 
	can be in state 0 (susceptible)	or in state 1 (infected).  In this model, 
	particles transition to infected at rate 	alpha.
	
	Inputs:
		alpha (float): Rate at which a particles react
		D (float): Diffusion coefficient
		L (float): Length of spatial domain
		S0 (int): Initial density of susceptible particles
		Nt (int): Number of timesteps
		dt (float): Length of timestep
		
	Uses:
		r_nhbd (float): radius of the neighborhood to look for infected particles
		dr (ndarray): Collect number of infected in a neighborhood of each susceptible
		id_xs (indices): Indices for which a particle is susceptible
		id_xi (indices): Indices for which a particle is infected
		xloc, yloc (float): location of a given susceptible particle
		p_decay (ndarray): Compute probability of decay
		R (ndarray): Compare to p_decay to determine whether S -> I
		ds (binary array): Change states if S -> I

	Outputs:
		t (ndarray): Times
		x (ndarray): x-Trajectories of Np particles
		y (ndarray): y-Trajectories of Np particles
		s (ndarray): States of Np particles
	"""
	# -- Number of particles is the density times the area -
	Np = L**2*S0				
	
	# -- Initialize time, position, state arrays -
	t = np.zeros( Nt+1 )
	x = np.zeros( (Nt+1, Np) )
	y = np.zeros( (Nt+1, Np) )
	s = np.zeros( (Nt+1, Np) )

	# -- Set random positions with one infected at center -
	x[0, :] = rng.uniform(0, L, Np)
	y[0, :] = rng.uniform(0, L, Np)
	x[0, 0] = L/2
	y[0, 0] = L/2
	s[0, 0] = 1
	
	# -- Go through time steps -
	for kt in range(Nt):
		# -- Detect infected particles (state 1) within x_nbhd of susceptible (state 0) -
		r_nbhd = .1
		dr = np.zeros(Np)
		id_xs = np.nonzero(1-s[kt,:])
		id_xi = np.nonzero(s[kt,:])
		for k_xs in id_xs[0]:
			xloc = x[kt, k_xs]
			yloc = y[kt, k_xs]
			r_check = (xloc - x[kt, id_xi])**2 + (yloc - y[kt, id_xi])**2
			dr[k_xs] = np.count_nonzero( r_check < r_nbhd**2)
		# -- Note that we need to divide by an area to ensure a density (units matter!) -
		p_decay = dr*dt*alpha/(np.pi*r_nbhd**2)

		# -- Determine which states switch -
		R = rng.uniform(0, 1, Np)
		ds = (p_decay > R)*1
		
		# -- Update time, state vector -
		s[kt+1, :] = s[kt, :] + ds
		t[kt+1] = t[kt] + dt

		# -- Move particles -
		dx = np.sqrt(2*D*dt)*rng.normal(0, 1, Np)
		dy = np.sqrt(2*D*dt)*rng.normal(0, 1, Np)
		x[kt+1,:] = x[kt, :] + dx
		y[kt+1,:] = y[kt, :] + dy

		# -- Ensure reflecting Boundary
		x[kt+1,:] = np.abs(x[kt+1,:])
		y[kt+1,:] = np.abs(y[kt+1,:])
		x[kt+1,:] = L - np.abs(L - x[kt+1,:])
		y[kt+1,:] = L - np.abs(L - y[kt+1,:])
	return s, t, x, y

# =============================================================================
# Main Simulation Function
# =============================================================================
def agent_SIR_2d():
	"""
	Identify the parameters needed, run the simulation, do some plotting
	"""

	# -- Parameters -
	alpha = 1	# Rate constant for S -> I
	D = 1		# Diffusion coefficient
	L = 10		# Length of spatial interval
	S0 = 1		# Density of S (particles per length) on 0 < x < L
	Nt = 3000	# Number of time steps
	dt = 0.01	# Length of time step

	s, t, x, y = xyReactDiff(alpha, D, L, S0, Nt, dt)

	# -- Plot trajectories as a movie-
	fig, ax = plt.subplots()
	xs_data, ys_data = [], []
	xi_data, yi_data = [], []
	s_plt = ax.plot([], [], 'g.')
	i_plt = ax.plot([], [], 'r.')
	
	def init():
		ax.set_xlim(0,L)
		ax.set_ylim(0,L)
		return s_plt, i_plt,

	def update(frame):
		ks = np.nonzero(1 - s[:, frame])
		ki = np.nonzero(s[:, frame])
		ps_update.set_xdata(x[ks, frame])
		ps_update.set_ydata(y[ks, frame])
		pi_update.set_xdata(x[ki, frame])
		ps_update.set_xdata(y[ki, frame])
		ax.set(title=f'Time t = {t[frame]:.3f} s')
		return ps_update, pi_update
		
		

	"""	
	plt.figure()
	frame = 13
	ks = np.nonzero(1 - s[:, frame])
	ki = np.nonzero(s[:, frame])
	plt.plot(x[ks, frame], y[ks, frame], 'g')
	plt.plot(x[ki, frame], y[ki, frame], 'r')
	plt.show()
	"""

	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt-1), interval=10)
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    agent_SIR_2d()

