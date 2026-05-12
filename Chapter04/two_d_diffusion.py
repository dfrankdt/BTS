#!/usr/bin/env python3
"""
Two Dimensional Diffusion

We extend the stochastic path from single_particle_diffusion.py to two dimensions.
The path simulates Brownian motion by solving a stochastic differential equation.

Figures produced
 - Figure 1: Final location of Np particles over Nt steps with step size dt
 - Figure 2: Probability of particle end position within a circle of radius R
 - Figure 3: Comparison of the probability distribution
 - Figure 4: Mean-squared displacement as a function of time
 
We note that the figures mirror Figures 4.5 and 4.6 in the text. Like the text,
we do not address how those results change as a function of D, leaving that 
instead for the exercises.

Just adding this line to see
 
This script is based on two_d_diffusion.m
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()

# =============================================================================
# Stochastic Differential Equation
# =============================================================================
def xyStochastic(D, dt, Nt, Np):
	"""
	Step forward with Np particles

		X_t+1 = X_t + dX_t
      
	where dX_t = sqrt(2*D*dt) N(0, 1, (Np, 2)) for some small dt.
    
	Inputs:
		D (float): Diffusion coefficient
		dt (float): Temporal step
		Nt (int): Number of steps
		Np (int): Number of particles
        
	Uses:
		kt (int): Temporal counter
		dX (ndarray): Spatial step
		r (float): Distance from the origin
		Nt (int): Large number to prepopulate
        
	Returns:
		t (ndarray): Time steps
		X (ndarray): Trajectory of Np particles
	"""


	X = np.zeros( (Np, 2, Nt+1) )
	t = np.zeros(Nt+1)
	for kt in range(Nt):
		dX = np.sqrt(2*D*dt)*rng.normal(0, 1, (Np, 2)) 
		X[:, :, kt+1] = X[:, :, kt] + dX
		t[kt + 1] = t[kt] + dt

	return t, X

# =============================================================================
# Main Simulation
# =============================================================================
def two_d_diffusion():

	# --- Global Parameters
	tf = 10			# End time
	dt = 1/2**10	# Step size
	Nt = int(tf/dt)	# Number of steps

	Np = 1000		# Number of particles
	D = 1			# Diffusion Coefficient

	# --- Collect and plot the end positions (Figure 1)
	t, X = xyStochastic(D, dt, Nt, Np)
	x, y = X[:,0,:], X[:,1,:]
	fig1, ax1 = plt.subplots()
	ax1.plot(x[:,-1], y[:,-1],'.b')
	ax1.set(xlabel = 'x', ylabel = 'y')
	ax1.set(title = rf'End Position of $N = {Np:1d}$ Two Dimensional Random Walks')

	# --- Compute the distribution outside a circle of radius r
	r0 = 15
	Nr = 2**5
	r = np.linspace(0, r0, Nr + 1)
	P = np.zeros(np.size(r))

	for kr in range(Nr+1):
		for kp in range(Np):
			P[kr] += (x[kp,-1]**2 + y[kp,-1]**2 < r[kr]**2)
	P = P/Np			

	# --- Plot the distribution P(R) (Figure 2)
	fig2, ax2 = plt.subplots()
	ax2.plot(r**2, P)
	ax2.set(xlabel = r'$R^2$', ylabel = r'$P(R)$')
	ax2.set(title = r'Fraction of points $r < R$')

	# --- Compare the log of the distribution to -R^2/40 (Figure 3)
	fig3, ax3 = plt.subplots()
	ax3.plot(r**2, -1/40*r**2, '--r', label = 'Theoretical Result')
	ax3.plot(r**2, np.log(1-P), 'k', label = 'Stochastic Result')
	ax3.set(xlabel = r'$R^2$', ylabel = r'$\log(1 - R)$')
	ax3.legend(loc = 'upper right')

	# --- Compute and plot the mean squared displacement
	R2 = x**2 + y**2
	R2ms = np.mean(R2, axis=0)
	fig4, ax4 = plt.subplots()
	ax4.plot(t, 4*t, '--r', label='Theoretical')
	ax4.plot(t, R2ms, '-k', label='Actual')
	ax4.set(xlabel = 't', ylabel = 'Mean Square Displacement')
	ax4.set(title = rf'Result of $N = {Np:1d}$ Two-Dimensional Random Walks')
	ax4.legend(loc='upper left')
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    two_d_diffusion()


