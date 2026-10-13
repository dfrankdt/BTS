#!/usr/bin/env python3
"""
CN with Upwind for Amobae and cAMP: We use an upwind scheme on the flux term
and CN on the second derivative to find a numerical approximation to the 
solution of the RD system.

TO DO: Something is wrong with my upwinding. Ugh.
"""

# =============================================================================
# Packpages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation
rng = np.random.default_rng()

# =============================================================================
# Numerical Solve
# =============================================================================
def doSolve(x, Du, Dc, rho, K, chi, uinit, cinit):
	# --- Spatial Discretization
	dx = x[1] - x[0]
	Nx = len(x) - 1
        
	# --- Spatial midpoints
	xm = np.zeros(Nx+2)
	xm[1:Nx+1] = (x[0:Nx] + x[1:Nx+1])/2
	xm[0] = x[0] - dx
	xm[-1] = x[-1] + dx
        
	# --- Temporal Discretization
	Nt = 2**5
	Nt_steps = 2**11
	dt = dx
	t = np.zeros(Nt+1)
        
	# --- State variables
	U = np.zeros( (Nx+1, Nt+1) )
	U[:, 0] = uinit
	C = np.zeros( (Nx+1, Nt+1) )
	C[:, 0] = cinit
	
	# --- CN Structure
	gamu = Du*dt/(2*dx**2)
	gamc = Dc*dt/(2*dx**2)
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	# --- No Flux
	D2[0, 0], D2[-1, -1] = -1, -1
	
	Acnu = np.eye(Nx+1) - gamu*D2
	Acnc = np.eye(Nx+1) - gamc*D2
	Bcnu = np.eye(Nx+1) + gamu*D2
	Bcnc = np.eye(Nx+1) + gamc*D2
	
	# --- Double sum to include intermediate steps
	uk, ck = uinit, cinit
	tk = 0
	for kt in range(Nt):
		for kt_inter in range(Nt_steps):
			# --- Flux for upwinding
			cplus = np.zeros(Nx+2)
			cplus[:-1] = ck
			cminus = np.zeros(Nx+2)
			cminus[1:] = ck
			chalf = (cplus + cminus)/2
			v = chi * (chalf[1:] - chalf[:-1])/dx
			Jplus = np.zeros(Nx+2)
			Jplus[1:] = v * uk * (v>0)
			Jminus = np.zeros(Nx+2)
			Jminus[:-1] = v * uk*(1 - (v>0))
			J = ( Jplus + Jminus )
			dJ = (J[1:] - J[:-1])/dx
			
			# --- Do the CN
			yu = Bcnu@uk - dt * dJ
			yc = Bcnc@ck + dt * (rho*uk - K*ck)
			ukp1 = np.linalg.solve(Acnu, yu)
			ckp1 = np.linalg.solve(Acnc, yc)
			
			uk = ukp1
			ck = ckp1
			tk += dt
		U[:, kt+1] = ukp1
		C[:, kt+1] = ckp1
		t[kt+1] = tk
	return t, U, C


# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_upwind_chemotax():
	# --- Parameters
	L = 10 		# length of the domain
	Du = 0.07 	# diffusion coefficient for dd
	Dc = 0.25 	# diffusion coefficient of cAMP
	rho = 0.2 	# rate of cAMP production
	K = 0.3 	# rate of cAMP degradation
	chi = 0.2 	# chemotactic parameter

	# --- Discretization
	Nx = 2**8
	x = np.linspace(0, L, Nx + 1)
	
	# --- Initialization
	u0 = np.ones(Nx+1) + 0.1/np.cosh(10 * (x - L/3.9) )
	c0 = rho/K*u0
	
	t, U, C = doSolve(x, Du, Dc, rho, K, chi, u0, c0)
	
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.8))
	ax1.plot(x, u0, '--')
	ax2.plot(x, c0, '--')
	ax1.plot(x, U[:, 1:])
	ax2.plot(x, C[:, 1:])
	ax1.set(xlabel = 'x', ylabel = 'u(x, t)', ylim = (0, 18))
	ax2.set(xlabel = 'x', ylabel = 'c(x, t)', ylim = (0, 18))
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_upwind_chemotax()
