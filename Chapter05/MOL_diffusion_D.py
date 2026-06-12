#!/usr/bin/env python3
"""
Diffusion with Dirichlet BCs via the Method of Lines

We approximate the solution of the Diffusion Equation (5.37) -- (5.38) with
Dirichlet boundary conditions using the method of lines. We discretize the
spatial region into N subintervals.  For the  N-1 interior points, we create
the second difference operator. For the two boundary points, we set the RHS
of the ODE to be zero.

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
def de_rhs(t, u, alpha):
	N = len(u)-1
	d2u = np.zeros(N+1)
	d2u[1:N] = u[0:N-1] - 2*u[1:N] + u[2:N+1]
	d2u[0] = 0
	d2u[-1] = 0
	du = alpha*d2u
	return du

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
	u0, uL = 0, 0	# Dirichlet Conditions

	# --- Discretization
	Nt, Nx = 2**8, 2**6
	dx = L/Nx
	dt = 0.001
	tf = dt*Nt
	x = np.linspace(0, L, Nx + 1)
	t = np.linspace(0, tf, Nt + 1)

	# --- Initial profile (Gaussian)
	uinit = np.exp(- (x - L/2)**2/(2*dx) )
	uinit[0] = u0
	uinit[-1] = uL

	# --- Initialization
	alpha = D/dx**2
	U = np.zeros( (Nx+1, Nt+1) )

	# --- Solve the IVP
	soln = solve_ivp(de_rhs, [0, tf], uinit, args=[alpha], dense_output=True)
	U = soln.sol(t)

	# --- Animation
	doMovie(x, t, U, 2**3)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	MOL_diffusion_D()
