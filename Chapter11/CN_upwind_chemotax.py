#!/usr/bin/env python3
"""
CN with Upwind for Amobae and cAMP: We use an upwind scheme on the flux term
and CN on the second derivative to find a numerical approximation to the 
solution of the RD system.

Produces
 - Figure 1: Snapshots of chemotactic amobae and chemoattractant, see Figure 11.3
 - Figure 2: Animation of snapshots shown in Figure 1
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# Numerical Solve
# =============================================================================
def doSolve(x, Du, Dc, rho, K, chi, uinit, cinit):
	# --- Spatial Discretization
	dx = x[1] - x[0]
	Nx = len(x) - 1
        
	# --- Temporal Discretization
	Nt = 2**4			# for outer loop
	t = np.zeros(Nt+1)	# times to record solution
	dt = dx				# step size for recording
	Nt_steps = 2**11	# intermediate steps
        
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
			# --- Half steps for c
			chalf = (ck[1:] + ck[:-1])/2
			
			# --- Velocity and upwind flux
			v = chi * (ck[1:] - ck[:-1])/dx
			J = np.zeros(Nx+2)
			J[1:-1] = v * (uk[:-1] * (v>0) + uk[1:] * (1 - (v>0)))
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
# Do the plotting
# =============================================================================
def doPlots(x, U, C):
	# --- Initialization
	u0, c0 = U[:, 0], C[:, 0]
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.8))
	ax1.plot(x, u0, '--')
	ax2.plot(x, c0, '--')
	
	# --- Plot the following timesteps
	ax1.plot(x, U[:, 1:])
	ax2.plot(x, C[:, 1:])
	ax1.set(xlabel = 'x', ylabel = 'u(x, t)', ylim = (0, 18))
	ax2.set(xlabel = 'x', ylabel = 'c(x, t)', ylim = (0, 18))
	return fig, ax1, ax2

# =============================================================================
# Animation
# =============================================================================
def doMovie(x, t, U, C):
	# --- Initialization
	Nt = len(t) - 1
	u0, c0 = U[:, 0], C[:, 0]
	
	# --- Initialize movie
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.8))
	u_plot_init = ax1.plot(x, u0, '--r', label = 'Initial Profile')
	c_plot_init = ax2.plot(x, c0, '--r', label = 'Initial Profile')
	u_update = ax1.plot([], [], '-b', label = 'Time Evolution')[0]
	c_update = ax2.plot([], [], '-b', label = 'Time Evolution')[0]
	ax1.set(xlabel = 'x', ylabel = 'u(x, t)', ylim=(0, 18))
	ax2.set(xlabel = 'x', ylabel = 'c(x, t)', ylim=(0, 18))
	ax1.legend(loc = 'upper right')
	ax2.legend(loc = 'upper right')
	
	# --- Funciton to update plot
	def update(frame):
		tk = t[frame]
		uk = U[:, frame]
		ck = C[:, frame]
		u_update.set_xdata(x)
		u_update.set_ydata(uk)
		c_update.set_xdata(x)
		c_update.set_ydata(ck)
		ax1.set(title=f'Chemotatic Amobae, t = {tk:4.1f} s')
		ax2.set(title=f'Chemoattractant, t = {tk:4.1f} s')
		return (u_update, c_update)
	
	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt+1), interval=100)
	return ani

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
	
	# --- Initialization and solve
	u0 = np.ones(Nx+1) + 0.1/np.cosh(10 * (x - L/3.9) )
	c0 = rho/K*u0
	t, U, C = doSolve(x, Du, Dc, rho, K, chi, u0, c0)
	
	# --- Plotting and Animation
	fig, ax1, ax2 = doPlots(x, U, C)
#	plt.show()

	ani = doMovie(x, t, U, C)
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_upwind_chemotax()
