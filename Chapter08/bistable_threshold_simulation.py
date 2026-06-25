#!/usr/bin/env python3
"""
Bistable Threshold Simulation

We identify the threshold behavior in the ODE associated with the traveling 
wave. We use Crank-Nicolson to solve the PDE and plot frames based on certain.

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# =============================================================================
#  Nonlinearity
# =============================================================================
def F(u, alpha):
	y = u * (1 - u) * (u - alpha)
	return y

# =============================================================================
# Crank-Nicolson Method
# =============================================================================
def doCN(x, t, uinit, alpha):
	""" 
	Crank-Nicolson to simulate the Bistable Equation
		
		u_t = u_xx +  u (1 - u) (u-a)
		
	with no-flux boundary conditions.
	"""
	dx = x[1]-x[0]
	dt = t[1]-t[0]
	Nx = np.size(x) - 1
	Nt = np.size(t) - 1

	U = np.zeros( (Nx+1, Nt+1) )
	U[:,0] = uinit

	# -- Matrices for performing CN
	gam = dt/(2*dx**2)
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	# -- No Flux BCs
	D2[0,0], D2[-1, -1] = -1, -1

	Acn = np.eye(Nx+1) - gam*D2
	Bcn = np.eye(Nx+1) + gam*D2
	
	# -- Initialization of CN method
	uk = uinit
	for kt in range(Nt):
		y = Bcn@uk + dt*F(uk, alpha)
		ukp1 = np.linalg.solve(Acn, y)
		U[:,kt+1] = ukp1
		uk = ukp1
	return U

# =============================================================================
# Create plots
# =============================================================================


# =============================================================================
# Main Simulation Function
# =============================================================================
def bistable_threshold_simulation():
	# --- Parameters
	L = 30			# Spatial Domain
	Tf = 100		# End time
	alpha = 0.25 	# Third zero of the cubic

	# --- Spatial and Temporal Scales 
	Nt, Nx = 2**10, 2**7
	dt, dx = Tf/Nt,  L/Nx
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, dt*Nt, Nt+1)

	# --- no traveling wave
	a = 0.40
	lam = 4.2
	u0_profile = a*(1/np.cosh(x/lam))**2
	U = doCN(x, t, u0_profile, alpha)

	fig, ax = plt.subplots()
	ax.plot(x, u0_profile, '--r')
	ax.set(xlabel = r'$\xi$', ylabel = r'u($\tau$, $\xi$)')
	ax.set(ylim = (-0.1, 1.1))
	for kt in range(0, Nt+1, 2**6):
		ax.plot(x, U[:, kt])
	plt.show()

	# --- yes traveling wave
	a = 0.41
	lam = 4.2
	u0_profile = a*(1/np.cosh(x/lam))**2
	U = doCN(x, t, u0_profile, alpha)

	fig, ax = plt.subplots()
	ax.plot(x, u0_profile, '--r')
	ax.set(xlabel = r'$\xi$', ylabel = r'u($\tau$, $\xi$)')
	ax.set(ylim = (-0.1, 1.1))
	for kt in range(0, Nt+1, 2**6):
		ax.plot(x, U[:, kt])
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    bistable_threshold_simulation()
