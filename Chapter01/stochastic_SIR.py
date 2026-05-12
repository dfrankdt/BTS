#!/usr/bin/env python3
"""
SIR Model via Gillespie Algorithm 

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
rng = np.random.default_rng()

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
			R = rng.random(2)
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
  	"""
	# --- Global Parameters
	alpha = 1
	N = 50
	Ntrials = 2000

	# -- Experiment 1: R0 = 2.5
	beta = 20
	t1 = np.zeros(Ntrials)
	ns1 = np.zeros(Ntrials)
	for k in np.arange(Ntrials):
		ns1[k], t1[k] = ustochastic(alpha, beta, N)

	figa, axa = plt.subplots()
	axa.plot(ns1, t1, '.')
	axa.set(xlabel = 'Number of Survivors', ylabel = 'Recovery Time')
	axa.set(title = r'Survivors and Recovery Time for $R_0 = 2.5$')

	# -- Experiment 2: R0 = 0.9
	beta = 55
	t2 = np.zeros(Ntrials)
	ns2 = np.zeros(Ntrials)
	for k in np.arange(Ntrials):
		ns2[k], t2[k] = ustochastic(alpha, beta, N)

	figb, axb = plt.subplots()
	axb.plot(ns2, t2, '.')
	axb.set(xlabel = 'Number of Survivors', ylabel = 'Recovery Time')
	axb.set(title = r'Survivors and Recovery Time for $R_0 = 0.9$')

	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    stochastic_SIR()

