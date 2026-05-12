#!/usr/bin/env python3

# --- TO DO: Book uses $\alpha t$ on the horizontal axis

"""
Decay Process via Gillespie Algorithm and ODE Simulation

We establish a decay process via the Gillespie Algorithm and compare the result
with the solution of an ODE system describing the evolution of the probability
distribution.

Reaction:
	S -> I with rate alpha
	
Figures Produced
 - Figure 1: 
 - Figure 2:
 
Note that this code is based on exponential_decay_via_Gillespie.m
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
rng = np.random.default_rng()

# =============================================================================
# Deterministic solution
# =============================================================================
def usoln(t, alpha, s0):
	"""
    Compute the solution of the related deterministic system
    
      s' = -alpha s
      s(0) = s0

    Inputs:
	    t (ndarray): Time
        alpha (float): Decay rate
        s0 (float): Initial number of particles
        
    Returns:
        s (ndarray): Array of solution
    """
	s = s0 * np.exp(-alpha * t)
	return s

# =============================================================================
# Stochastic solution
# =============================================================================
def ustochastic(alpha, N):
	"""
	Perform the Gillespie algorithm, computing the times t at which the transitions
	from k to k - 1 particles occur.

	Inputs:
		alpha (float): Decay rate
		N (int): Initial number of particles

	Returns:
		t (ndarray): times at which a transitions occur
		k (int array): number of particles at the time t
	"""
	R = rng.random(N - 1)
	k = np.arange(N, 0, -1)
	t = np.zeros(N)
	t[1:] = -1/alpha * np.cumsum( np.log(R) / (k[1:] ))
	return t, k

# =============================================================================
# RHS of System of Stochastic DEs
# =============================================================================
def de_rhs(t, p, alpha):
    """
    Compute the right hand side of the differential equation for p_k(t), given by
        pN' = - alpha * (N) * pN
        pk' = - alpha * (k) * pk + alpha * (k+1) * pk+1
        p1' =   alpha * (2) * p2
    for k = 1, 2, ..., N-1 

    Inputs:
        t (float): Time (unused since the DE is autonomous)
        p (ndarray): Array of probabilities for k = 0, 1, ..., N 
        alpha (float): Decay rate

    Returns:
        dp (ndarray): Array of time derivatives for each pk
    """
    N = np.size(p) 
    k = np.arange(0, N)
    dp = np.zeros(N)
    dp[0] = alpha * k[1] * p[1]
    dp[1:-1] = - alpha * (k[1:-1]) * p[1:-1] + alpha * (k[2:]) * (p[2:])
    dp[-1] = - alpha * (k[-1]) * p[-1]
    return dp

# =============================================================================
# Main Simulation Function
# =============================================================================
def exponential_decay_via_Gillespie():
	N = 25			# Particles
	alpha = 1		# Decay/death rate
	Ntrials = 10000	# Independent trials

    # --- One stochastic trajectory via Gillespie 
	[tzero, Ndecay] = ustochastic(alpha, N)
	
	# --- Solution of deterministic system
	tf = tzero[-1]
	td = np.linspace(0, tf, 2**9+1)
	ud = usoln(td, alpha, N)
	
	# --- Create subfigure (a)
	figa, axa = plt.subplots()
	axa.plot(td, ud, '--r')
	axa.step(tzero, Ndecay, where='post')
	axa.set(xlabel = 't', ylabel = 'Number of Particles')
	axa.set(title = 'Decay Trajectory: Stochastic vs. Deterministic')

    # --- Independent Gillespie trials
	tlist = np.zeros( (Ntrials, N) )
	for k in range(Ntrials):
		[tlist[k,], Ndecay] = ustochastic(alpha, N)
	t0dist, bins = np.histogram(tlist[:,-1], 50, density=True)
	bin_center = (bins[:-1] + bins[1:])/2
	
	# --- Compare the actual and predicted mean extinction time
	tmean = np.mean(tlist[:,-1])
	tmean_pred = 1/(alpha)*(np.log(N) + np.euler_gamma)
	print('The mean time until total extinction', tmean)
	print('The predicted mean time until extinction is', tmean_pred)
        
	# --- Solution of the ODE system for pk
	tf = max(tlist[:,-1])
	pinit = np.zeros(N+1)
	pinit[-1] = 1

	soln = solve_ivp(de_rhs, [0, tf], pinit, args=[alpha], dense_output=True)
	tt = np.linspace(0, tf, 2**9+1)
	p = soln.sol(tt).T
	dp0 = alpha * p[:,1]

	# --- Create subfigure (b)
	figb, axb = plt.subplots()
	axb.plot(tt, dp0, '--r')
	axb.plot(bin_center, t0dist, '.')
#	axb.stairs(t0dist, bins, fill=True)
	axb.set(xlabel = 'Extinction Time', ylabel = r'$dp_0/dt$')
	axb.set(title ='Total Decay: Data vs ODE Prediction')
	
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    exponential_decay_via_Gillespie()
