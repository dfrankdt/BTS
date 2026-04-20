#!/usr/bin/env python3
"""
SIR Model via Gillespie Algorithm and ODE Simulation

This script simulates the SIR Model through a stochastic process
(simulated with the Gillespie algorithm) and compares results for
cases R0 < 1 and R0 > 1

Reaction:
    S + I -> 2I  with rate alpha
	I -> R with rate beta

Generates figures:
  - Figure 1: Time to clear infection and number of infected for R0 = 2.5
  - Figure 2: Time to clear infection and number of infected for R0 = 0.9

Note that this code is based on stochastic_SIR.m
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# Stochastic solution
# =============================================================================
def ustochastic(alpha, beta, N):
	"""
	Perform the Gillespie algorithm, computing the times t at which the transitions
	between states occur.
	        
	Inputs:
		alpha (float): Infection rate
		beta (float): Recovery rate
		N (int): Initial number of particles
	
	Uses:
		u (ndarray): [ns, ni, nr], numbers of SIR
		C (ndarray): Change matrix 
		t (float): time
		R (ndarray): random numbers to determine wait time
		h (ndarray): array containing reaction rates
		dt (float): amount of time to the next reaction
		z (ndarray): array to relay which reaction occurs
                
	Returns:
        t (ndarray): times at which a transitions occur
		k (int array): number of particles at the time t
	"""
	u = np.array([N, 1, 0])  
	C = np.array([ [-1, 1, 0], [0, -1, 1] ])
	t = 0
	# -- Step forward while there is a nonzero infected population
	while u[1] > 0:
		if u[0] > 0:
			R = np.random.rand(2)
			h = np.array([alpha*u[0]*u[1], beta*u[1]])
			dt = min(-np.log(R)/h)
			z = np.where(-np.log(R)/h > dt, np.zeros(np.size(R)), np.ones(np.size(R)))
		# -- If no susceptibles just step forward according to infected reaction
		else:
			dt = -np.log(R[1])/h[1]
			z = np.array([0, 1])
		t = t + dt
		u = u + z@C
	return u[0], t
    	

# =============================================================================
# Main Simulation Function
# =============================================================================
def stochastic_SIR():
	"""
	Simulate stochastic SIR model via the Gillespie algorithm with values alpha
	and beta designed to produce relative reproductive number R0 = 2.5 and R0 = 0.9

	Generates figures:
	  - Figure 1: Time to clear infection and number of infected for R0 = 2.5
	  - Figure 2: Time to clear infection and number of infected for R0 = 0.9
  	"""

	# -- Experiment 1: R0 = 2.5
	kfig = 1
	alpha = 1
	beta = 20
	N = 50
	ktrials = 2000
	t = np.zeros(ktrials)
	ns = np.zeros(ktrials)
	for k in np.arange(ktrials):
		ns[k], t[k] = ustochastic(alpha, beta, N)

	plt.figure(kfig)
	plt.plot(ns, t, '.')
	plt.xlabel('Number of Survivors')
	plt.ylabel('Recovery Time')

	# -- Experiment 2: R0 = 0.9
	kfig = 2
	alpha = 1
	beta = 55
	N = 50
	ktrials = 2000
	t = np.zeros(ktrials)
	ns = np.zeros(ktrials)
	for k in np.arange(ktrials):
		ns[k], t[k] = ustochastic(alpha, beta, N)

	plt.figure(kfig)
	plt.plot(ns, t, '.')
	plt.xlabel('Number of Survivors')
	plt.ylabel('Recovery Time')

	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    stochastic_SIR()

