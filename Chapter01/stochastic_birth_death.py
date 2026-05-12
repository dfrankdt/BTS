#!/usr/bin/env python3
"""
Stochastic Birth and Death via Gillespie Algorithm and ODE Simulation

We establish a birth/death process via the Gillespie Algorithm and compare the 
result with the solution of an ODE system describing the evolution of the 
probability distribution.

Reaction:
    S -> I  with rate alpha
    S -> 2S with rate beta

Figures Produced
 - Figure 1: 
 - Figure 2:
 
Note that this code is based on stochastic_birth_death.m
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
rng = np.random.default_rng()

# =============================================================================
# Deterministic Solution
# =============================================================================
def usoln(t, alpha, beta, s0):
	"""
	The solution of the deterministic differential equation (when alpha > beta)
	is a decaying exponential.
	"""
	r = alpha - beta
	s = s0*np.exp(-r*t)
	return s

# =============================================================================
# Stochastic Solution
# =============================================================================
def ustochastic(alpha, beta, N):
	"""
	Perform the Gillespie algorithm, computing the times t at which the transitions
	between states occur.

	Inputs:
		alpha (float): Death rate
		beta (float): Birth rate
		N (int): Initial number of particles

	Uses:
		u (ndarray): Number of S (artifically populated to size 10000)
		t (float): time
		lam (ndarray): Transition rates between states
		R (ndarray): random numbers to determine wait time and which reaction
		h (ndarray): array containing reaction rates
		dt (float): amount of time to the next reaction
		z (ndarray): array to relay which reaction occurs
	"""
	t = np.zeros(1000)
	u = np.zeros(1000)
	u[0] = N

	kt = 0
	while u[kt]>0:
		lam = np.array([u[kt]*alpha, u[kt]*beta])
		r = sum(lam)
		x = np.cumsum(lam)/r
		R = rng.random(2)
		
		dt = -np.log(R[0])/r
		du = (-1)*(R[1]<x[0]) + (1)*(R[1]>x[0])
		t[kt+1] = t[kt] + dt
		u[kt+1] = u[kt] + du
		kt += 1
	return t[:kt+1], u[:kt+1]

def de_rhs(t, p, alpha, beta):
    """
    Compute the right hand side of the differential equation for p_k(t), given by
        pN' = - alpha * (N) * pN + beta * (N-1) * pN-1
        pk' = - (alpha * (k) + beta *(k)) * pk + alpha * (k+1) * pk+1 + beta * (k-1) * pk-1
        p0' =   alpha * (1) * p1
    for k = 1, 2, ..., N-1. Note that we artifically truncate the system at N particles.

    Inputs:
        t (float): Time (unused since the DE is autonomous)
        p (ndarray): Array of probabilities for k = 0, 1, ..., N 
        alpha (float): Death rate
        beta (float): Birth rate

    Returns:
        dp (ndarray): Array of time derivatives for each pk
    """
    N = np.size(p) 
    k = np.arange(N)
    dp = np.zeros(N)
    dp[0] = alpha * k[1] * p[1]
    dp[1:-1] = - (alpha * (k[1:-1]) + beta * (k[1:-1])) * p[1:-1]
    dp[1:-1] = dp[1:-1] + alpha * (k[2:]) * (p[2:]) + beta * k[:-2] * (p[:-2])
    dp[-1] = - alpha * (k[-1]) * p[-1] + beta * (k[-2]) * p[-2]
    return dp

# =============================================================================
# Main Simulation Function
# =============================================================================
def stochastic_birth_death():
	# --- Global Parameters
	N = 20					# Individuals
	alpha, beta = 1, 0.25	# Death and Birth rates
	Ntrials = 10000			# Independent trials

	# --- One stochastic trajectory via Gillespie
	t, u = ustochastic(alpha, beta, N)
	
	# --- Solution of deterministic system
	tf = t[-1]
	td = np.linspace(0, tf, 2**9+1)
	ud = usoln(td, alpha, beta, N)
	
	# --- Create subfigure (a)
	figa, axa = plt.subplots()
	axa.plot(td, ud, '--r')
	axa.step(t, u, where='post')
	axa.set(xlabel = 't', ylabel = 'Number of Individuals')
	axa.set(title = 'Birth and Death Trajectory: Stochastic vs Deterministic')

	# --- Independent Gillespie trials
	tlist = np.zeros(Ntrials)
	for k in range(Ntrials):
		t, u = ustochastic(alpha, beta, N)
		tlist[k] = t[-1]
	t0dist, bins = np.histogram(tlist, 50, density=True)
	bin_center = (bins[:-1] + bins[1:])/2
	tmean = np.mean(tlist)
	
	# --- Solution of ODE system for pk
	tf = max(tlist)
	pinit = np.zeros(5*N)
	pinit[N+1] = 1

	soln = solve_ivp(de_rhs, [0, tf], pinit, args=[alpha, beta], dense_output=True)
	tt = np.linspace(0, tf, 2**9+1)
	p = soln.sol(tt).T
	dp0 = alpha * p[:,1]
	
	# --- Create subfigure (b)
	figb, axb = plt.subplots()
	axb.plot(tt, dp0, '--r')
	axb.plot(bin_center, t0dist, '.')
#	axb.stairs(t0dist, bins, fill=True)
	axb.set(xlabel='Extinction Time', ylabel=r'$dp_0/dt$')
	axb.set(title='Birth and Death: Data vs ODE Prediction')

	plt.show()
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    stochastic_birth_death()
