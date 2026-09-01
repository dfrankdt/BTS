#!/usr/bin/env python3
"""
Agent based 2-dimensional SIR

This script simulates agent based process in which particles move randomly
but may interact according to S -> I with rate alpha if the infected particle
is close enough to a susceptible particle.

Produces an animation showing susceptible and infected

This script is based on agent_SIR_2d.py

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
def xyReactDiff(alpha, D, L, S0, r_nbhd, Nt, dt):
	"""
	Pass Np = (L**2 * S0) particles through Nt timesteps. Particles 
	can be in state 0 (susceptible)	or in state 1 (infected).  In this model, 
	particles transition to infected at rate alpha.
	
	Inputs:
		alpha (float): Rate at which a particles react
		D (float): Diffusion coefficient
		L (float): Length of spatial domain
		S0 (int): Initial density of susceptible particles
		r_nhbd (float): radius of the neighborhood to look for infected particles
		Nt (int): Number of timesteps
		dt (float): Length of timestep
		
	Uses:
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
	# --- Number of particles is the density times the area
	Np = L**2*S0
	
	# --- Initialize time, position, state arrays
	t = np.zeros( Nt+1 )
	x = np.zeros( (Nt+1, Np) )
	y = np.zeros( (Nt+1, Np) )
	s = np.zeros( (Nt+1, Np) )

	# --- Set random positions with one infected at center
	x[0, :] = rng.uniform(0, L, Np)
	y[0, :] = rng.uniform(0, L, Np)
	x[0, 0] = L/2
	y[0, 0] = L/2
	s[0, 0] = 1
	
	# --- Go through time steps
	for kt in range(Nt):
		# --- Detect infected particles (state 1) within r_nbhd of susceptible (state 0)
		dr = np.zeros(Np)
		id_xs = np.nonzero(1-s[kt,:])
		id_xi = np.nonzero(s[kt,:])
		for k_xs in id_xs[0]:
			xloc = x[kt, k_xs]
			yloc = y[kt, k_xs]
			r_check = np.sqrt((xloc - x[kt, id_xi])**2 + (yloc - y[kt, id_xi])**2)
			dr[k_xs] = np.count_nonzero( r_check < r_nbhd)
		# --- Note that we need to divide by an area to ensure a density (units matter!) 
		p_decay = dr*dt*alpha/(np.pi*r_nbhd**2)

		# --- Determine which states switch 
		R = rng.uniform(0, 1, Np)
		ds = (p_decay > R)*1
		
		# --- Update time, state vector
		s[kt+1, :] = s[kt, :] + ds
		t[kt+1] = t[kt] + dt

		# --- Move particles 
		dx = np.sqrt(2*D*dt)*rng.normal(0, 1, Np)
		dy = np.sqrt(2*D*dt)*rng.normal(0, 1, Np)
		x[kt+1,:] = x[kt, :] + dx
		y[kt+1,:] = y[kt, :] + dy

		# --- Ensure reflecting boundary
		x[kt+1,:] = np.abs(x[kt+1,:])
		y[kt+1,:] = np.abs(y[kt+1,:])
		x[kt+1,:] = L - np.abs(L - x[kt+1,:])
		y[kt+1,:] = L - np.abs(L - y[kt+1,:])
	return t, x, y, s

# =============================================================================
# Create Animation
# =============================================================================
def doMovie(t, x, y, s, L):
	# --- Initialize data structures and animation
	Nt = len(t) - 1
	dL = 0.1*L
	fig, ax = plt.subplots()
	s_plt = ax.plot([], [], '.g', label = 'Susceptible')[0]
	i_plt = ax.plot([], [], '.r', label = 'Infected')[0]
	ax.set(xlabel = 'x', ylabel = 'y')
	ax.set(xlim = (-dL, L+dL), ylim = (-dL, L+dL))
	ax.legend(handles = [s_plt, i_plt], loc='upper right')

	# --- Update in each frame
	def update(frame):
		tk = t[frame]
		# -- Get indices of susceptible and infected
		ks = np.nonzero(1 - s[frame, :])
		ki = np.nonzero(s[frame, :])
		Ni = np.count_nonzero(s[frame,:])
		# -- Update plots
		s_plt.set_xdata(x[frame, ks])
		s_plt.set_ydata(y[frame, ks])
		i_plt.set_xdata(x[frame, ki])
		i_plt.set_ydata(y[frame, ki])
		ax.set(title=f'Time t = {tk:.3f} s, Ni = {Ni:2d} Infected')
		return s_plt, i_plt

	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt+1), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def agent_SIR_2d():
	"""
	Identify the parameters needed, run the simulation, do some plotting
	"""

	# --- Parameters 
	alpha = 1	# Rate constant for S -> I
	r_nhbd = .1 # Radial neighborhood S -> I interaction occurs
	D = 1		# Diffusion coefficient
	L = 10		# Length of spatial interval
	S0 = 1		# Density of S (particles per length) on 0 < x < L
	Nt = 200	# Number of time steps
	dt = 0.01	# Length of time step

	t, x, y, s = xyReactDiff(alpha, D, L, S0, r_nhbd, Nt, dt)
	doMovie(t, x, y, s, L)

# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    agent_SIR_2d()

