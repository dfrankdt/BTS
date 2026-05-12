#!/usr/bin/env python3
"""
Discrete random walk

This script simulates a number of discrete random walks with probably of moving 2a

Figures produced
	- Figure 1: position after 100 time steps
	- Figure 2: mean squared displacement as a function of time step n compared to theory
	- Figure 3: comparison of histogram to Gaussian
	
That this code is based on discrete_random_walk.m
"""

# =============================================================================
# Packages
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()

# =============================================================================
# Main Simulation Function
# =============================================================================
def discrete_random_walk():
	"""
	Perhaps add some comments
	"""
	# --- Simulation parameters
	alpha = 0.1		# probability of moving left
	beta = alpha	# probability of moving right
	N = 2000		# trials
	ns = 100		# number of steps
	
	# --- Simulation tools
	X = np.zeros( (N, ns+1) )				# Positions
	tsteps = np.arange(ns+1)				# Steps 
	c = np.array([alpha, alpha + beta, 1]) 	# vector to determine direction
	xm = np.array([-1, 1, 0])				# vector to move L, R, to stay put
	
	# --- Pass through the ns steps
	for j in tsteps-1:
		R = rng.random(N)
		mL = np.where(R<c[0], np.ones(np.size(R)), np.zeros(np.size(R)))
		mR = np.where((R>c[0])*(R<c[1]), np.ones(np.size(R)), np.zeros(np.size(R)))
		mZ = np.where(R>c[1], np.ones(np.size(R)), np.zeros(np.size(R)))
		X[:,j+1] = X[:,j] + mL*xm[0] + mR*xm[1] + mZ*xm[2]

	# --- Compute the mean-squared distance, also the variance when alpha = beta
	Xmd2 = np.mean(X**2, 0)
	
	# --- Collect the end positions to plot against a Gaussian
	XeDist, bins = np.histogram(X[:,-1], bins='auto', density=True)
	bin_center = (bins[:-1] + bins[1:])/2
	z = np.linspace(-15, 15, 2**9+1)
	p = 1/(np.sqrt(4*np.pi*alpha*ns)) * np.exp(-z**2/(4*alpha*ns))
	
	# --- Do some plotting
	fig1, ax1 = plt.subplots()
	ax1.step(tsteps, X[:20,:].T)
	ax1.set(xlabel = 'Time Step', ylabel = 'Position')
	ax1.set(ylim=(-15, 15))
	
	fig2, ax2 = plt.subplots()
	ax2.plot(tsteps, 2*alpha*tsteps, '--r', label='Theoretical')
	ax2.plot(tsteps, Xmd2, '.k', label='Actual')
	ax2.set(xlabel = 'Time Step', ylabel = 'Mean Squared Displacement')
	ax2.legend(loc='upper left')
	
	fig3, ax3 = plt.subplots()
	ax3.plot(bin_center, XeDist, '.', label='Actual')
	#ax3.stairs(XeDist, bins, fill=True)
	ax3.plot(z, p, label='Theoretical')
	ax3.set(xlabel = 'End Position', ylabel = 'Density')
	ax3.legend(loc='upper right')
		
	plt.show()
		
# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    discrete_random_walk()


