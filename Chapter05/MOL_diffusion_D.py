#!/usr/bin/env python3
"""
Diffusion with Dirichlet BCs via the Method of Lines

We approximate the solution of the Diffusion Equation (5.37) -- (5.38) with
Dirichlet boundary conditions using the method of lines. We discretize the
spatial region into N subintervals.  For the  N-1 interior points, we create
the second difference operator. For the two boundary points, we set the RHS
of the ODE to be zero.

Specifically, we find a numerical solution to the BVP

 u_t = D u_xx
 u = U0 at x=0
 u = UL at x=L

where the nonhomogenous boundary conditions are given as Dirichlet conditions.
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
def de_rhs(t, u, dx, D):
	N = len(u)-1
	alpha = D/dx**2
	
	d2U = np.zeros(N+1)
	d2U[1:N] = u[0:N-1] - 2*u[1:N] + u[2:N+1]
	d2U[0] = 0
	d2U[-1] = 0
	dudt = alpha*d2U
	return dudt

# =============================================================================
# Animation
# =============================================================================
def doMovie(x, t, U, ktskip):
	# --- Initialize data structures
	uinit = U[:,0]
	Nt = len(t) - 1
		
	# --- Initialize movie
	fig, ax = plt.subplots()
	p_init = ax.plot(x, uinit, '--r', label='Initial Profile')
	p_update = ax.plot([], [], '-b', label='Time Evolution')[0]
	ax.set(xlabel='x', ylabel='u(x, t)')
	ax.legend(loc='upper right')

	# --- Function to update the plot with the current frame
	def update(frame):
		tk = t[frame]
		u = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(u)
		ax.set(title=f'Time t = {tk:.2f} s')
		return(p_update)

	ani = manimation.FuncAnimation(fig=fig, func=update, 
		frames=range(0, Nt+1, ktskip), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def MOL_diffusion_D():
	# --- Global parameters
	L = 1			# length of domain
	D = .1			# diffusion coefficient
	
	# --- Boundary Conditions
	U0, UL = 0, 0	# Dirichlet Conditions

	# --- Discretization
	Nx = 2**6		# number of spatial partitions
	Nt = 2**9		# number of temporal partitions
	dx = L/Nx		# spatial discretization
	dt = 0.001		# temporal discretization, dt < dx**2/(2D)
	tf = dt*Nt		# end time
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, tf, Nt+1) 

	# --- Initial profile, pass boundary conditions to solver
	IVP_arg = [dx, D]
	u0_profile = np.exp( -(x - L/2)**2/(2*dx) )
	u0_profile[0] = U0
	u0_profile[-1] = UL

	# --- Solve the IVP
	soln = solve_ivp(de_rhs, [0, tf], u0_profile, args=IVP_arg, dense_output=True)
	U = soln.sol(t)

	# --- Animation
	doMovie(x, t, U, 2**3)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	MOL_diffusion_D()
