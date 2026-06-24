#!/usr/bin/env python3
"""
Discrete diffusion via the Gillespie algorithm

This script uses a Gillespie algorithm to simulate diffusion of M particles 
in N boxes. We take N to be odd so that the entire initial distribution fits
in the middle box.

Figures produced:
 - Figure 1: 
 - Figure 2: 
 
This script is based on discrete_diffusion_via_Gillespie.m

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# =============================================================================
# Gillespie Algorithm
# =============================================================================
def uGillespie(M, N, Nt, alpha):
	"""
	Perform the Gillespie algorithm to simulate the diffusion of particles 
	between discrete boxes.
	
	Inputs:
		M (integer):	Total number of particles
		N (integer): 	Total number of boxes
		Nt (integer): 	Number of time steps
		alpha (float): 	Probability of leaving a box
	
	Returns:
		t (ndarray):	Array of time values
		U (ndarray):	Array of distribution over time steps
	"""
	U = np.zeros( (Nt+1, N) )
	t = np.zeros(Nt+1)

	# Create initial distribution 
	rng = np.random.default_rng()
	u = np.zeros(N)
	u[int((N-1)/2)] = M
	U[0,:] = u
	
	pa = alpha*np.ones(N)
	pa[1:-1] = pa[1:-1]*2

	for kt in range(Nt):
		R = rng.uniform(0, 1, 3)
		r = pa*u
		Rsig = np.sum(r)
		
		# --- Increment the time
		dt = -np.log(R[0])/Rsig

		# --- Do the moving
		mr = np.zeros(N)
		pj = r[0]/Rsig
		j = 0
		while pj<R[1]:
			j = j+1
			pj = pj + r[j]/Rsig
		# --- Ensure the particle in the first box moves right
		if j==0:
			mr[j], mr[j+1] = np.array([-1, 1])
		# --- Ensure the particle in the last box moves left
		elif j==N-1:
			mr[j-1], mr[j] = np.array([1, -1])
		# --- Choose right or left, depending 
		else:
			mr[j-1], mr[j], mr[j+1] = np.where(R[2]<0.5, np.array([0, -1, 1]), np.array([1, -1, 0]))
		u = u + mr
		U[kt+1,:] = u
		t[kt+1]  = t[kt] + dt
	return t, U

# =============================================================================
# Create animation 
# =============================================================================
def doMovie(x, t, U, M):
	u0 = U[0,:]
	Nt = len(t) - 1

	fig, ax = plt.subplots()
	p_update = ax.plot(x, u0/M, 'or', label='Distribution')[0]
	ax.set(xlabel='boxes', ylabel='proportion')
	ax.set(title='Distribution')
	ax.legend(loc='upper left')
	
	def update(frame):
		tk = t[frame]
		u = U[frame, :]/M
		p_update.set_ydata(u)
		ax.set(title=rf'Time $t = {tk:.3f}$ s')
		ax.legend(loc='upper left')
		return(p_update)

	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt+1), interval=10)
	plt.show()

# =============================================================================
# Main Simulation 
# =============================================================================
def discrete_diffusion_via_Gillespie():
	# --- Global Parameters
	M = 500 	# Particles
	N = 9   	# Boxes	
	alpha = 0.1	# Probability of moving L or R
	Nt = 1000	# Time steps

	tlist, U = uGillespie(M, N, Nt, alpha)

	# Create spatial scale to plot boxes
	x = np.arange(1, N+1)
	doMovie(x, tlist, U, M)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    discrete_diffusion_via_Gillespie()

