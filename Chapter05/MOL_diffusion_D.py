#!/usr/bin/env python3
"""
Diffusion with Dirichlet BCs via the Method of Lines

We approximate the solution of the Diffusion Equation (5.37) -- (5.38) with
Dirichlet boundary conditions using the method of lines. We discretize the
spatial region into N subintervals (N+1 endpoints) and consisder the N-1
interior points.

With zero Dirichlet conditions, the boundary points are fixed and we need only
solve the IVP on the interior points.

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
def de_rhs(t, u, alpha, u0, uL):
	N = len(u)-1
	d2u = np.zeros(N+1)
	d2u[1:N] = u[0:N-1] - 2*u[1:N] + u[2:N+1]
	d2u[0] = u0 - 2*u[0] + u[1]
	d2u[-1] = u[-2] - 2*u[-1] + uL
	du = alpha*d2u
	return du

# =============================================================================
# Animation
# =============================================================================
def doMovie(x, t, U):
	# --- Initialize data structures
	Nt = len(t) - 1
	uinit = U[:,0]
		
	# --- Initialize movie
	fig, ax = plt.subplots()
	p_init = ax.plot(x, uinit, '-r', label='Initial Profile')
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

	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt + 1), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def MOL_diffusion_D():
	# --- Global parameters
	D = 1			# diffusion coefficient
	L = 1			# length of domain
	Nx = 2**6		# number of partitions in x
	tf = 0.5		# end time
	dt = 0.01		# time step

	# --- Boundary Conditions
	u0, uL = 0, 0	        # Dirichlet Conditions

	# --- Variables
	x = np.linspace(0, L, Nx + 1)
	dx = L/Nx
	Nt = int(tf/dt)
	t = np.linspace(0, tf, Nt + 1)
	alpha = D/dx**2
	U = np.zeros( (Nx+1, Nt+1) )
	U[0, 1:] = u0*np.ones(Nt)
	U[-1, 1:] = uL*np.ones(Nt)

	# --- Initial profile (Gaussian)
	x_interior = x[1:Nx]
	uinit = np.exp(- (x_interior - L/2)**2/(2*dx) )

	# --- Solve the IVP on the interior points
	soln = solve_ivp(de_rhs, [0, tf], uinit, args=[alpha, u0, uL], dense_output=True)
	U[1:Nx,:] = soln.sol(t)

	# --- Animation
	doMovie(x, t, U)
    

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	MOL_diffusion_D()
