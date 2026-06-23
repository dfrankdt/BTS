#!/usr/bin/env python3
"""
Two-Dimensional Diffusion (radial symmetry) via MOL

We solve the two-dimensional diffusion equation, leveraging radial symmetry
using the Method of Lines.

TO DO: radially symmetric plots

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(t, y, R, D, U0):
	# --- Discretization structure on interior points
	N = len(y) + 1
	r = np.linspace(0, R, N+1)
	dr = R/N
	alpha = D/dr**2
	
	# --- Construct solution including boundary points
	u = np.zeros(N+1)
	u[1:N] = y
	u[0] = y[1]
	u[-1] = U0
	
	# --- Differences
	up = u[2:N+1] - u[1:N]
	um = u[1:N] - u[0:N-1]
	rp = (r[2:N+1] + r[1:N])/2
	rm = (r[1:N] + r[0:N-1])/2

	D2u = (up*rp - um*rm)/r[1:N]
	dudt = alpha*D2u
	return dudt

# =============================================================================
# Animation
# =============================================================================
def doMovie(r, t, U, ktskip):
	# --- Initialize data structures
	uinit = U[:,0]
	Nt = len(t) - 1
	Nr = len(r) - 1
	
	# --- Initialize movie
	fig, ax = plt.subplots()
	p_init = ax.plot(r, uinit, '--r', label='Initial Profile')
	p_update = ax.plot([], [], '-b', label='Time Evolution')[0]
	ax.set(xlabel='r', ylabel='u(r, t)')
	ax.legend(loc='upper right')
	
	# --- Function to update the plot with the current frame
	def update(frame):
		tk = t[frame]
		u = U[:, frame]
		p_update.set_xdata(r)
		p_update.set_ydata(u)
		ax.set(title=f'Time t = {tk:.2f} s')
		return(p_update)

	ani = manimation.FuncAnimation(fig=fig, func=update, 
		frames=range(0, Nt+1, ktskip), interval=100)

	plt.show()


	 

# =============================================================================
# Main Simulation Function
# =============================================================================
def Diffusion_2D_MOL():
	# --- Global Parameters
	R = 1	# Radius
	D = 1	# Diffusion cofficient
	
	# --- Boundary condition at r = 1
	U0 = 1
	
	# --- Discretization
	Nr = 2**6
	Nt = 2**6
	dr = R/Nr
	dt = 0.01
	tf = dt*Nt
	r = np.linspace(0, R, Nr+1)
	t = np.linspace(0, tf, Nt+1)
	
	# --- Initial profile, pass boundary conditions
	IVP_arg = [R, D, U0]
	u0_profile = np.zeros(Nr+1)
	u0_profile[-1] = U0
	
	# --- Solve the IVP on interior points
	U = np.zeros( (Nr+1, Nt+1) )
	U[:,0] = u0_profile
	soln = solve_ivp(de_rhs, [0, tf], u0_profile[1:Nr], args=IVP_arg, dense_output=True)
	U[1:Nr,:] = soln.sol(t)
	
	# --- Embed boundary conditions
	U[0,:] = U[1,:]
	U[-1, :] = U0
	
	# --- Animation?
	doMovie(r, t, U, 2**2)
	
	# --- Plotting
	fig1, ax1 = plt.subplots()
	u0 = U[0,:]
	ax1.plot(t, u0)
	ax1.set(xlabel = 't', ylabel = 'u(0, t)')
	
	fig2, ax2 = plt.subplots()
	lam = 5.77*D/R**2
	ax2.plot(t, np.log(1 - u0))
	ax2.plot(t, -lam*t, '--r')
	ax2.set(xlabel = 't', ylabel = r'$\ln$(1 - u(0, t))')
	
	plt.show()
	
	
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    Diffusion_2D_MOL()
