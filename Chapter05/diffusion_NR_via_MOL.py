#!/usr/bin/env python3
"""
Diffusion with Neumann or Robin BCs via the Method of Lines

We approximate the solution of the Diffusion Equation (5.37) -- (5.38) with
using either Neumann or Robin boundary conditions using the method of lines

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# DE RHS -- Neumann or Robin Condition
# =============================================================================
def de_rhs_NR(t, u, alpha, R, B):
	N = len(u)-1
	R0, RN = R
	B0, BN = B

	d2u = np.zeros(N+1)
	d2u[1:N] = u[0:N-1] - 2*u[1:N] + u[2:N+1]
	d2u[0] = (-2 - R0)*u[0] + 2*u[1] + B0
	d2u[-1] = 2*u[-2] + (-2 - RN)*u[-1] + BN
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
def diffusion_NR_via_MOL():
	# --- Global parameters
	D = 1			# diffusion coefficient
	L = 1			# length of domain
	Nx = 2**6		# number of partitions in x
	tf = 0.5		# end time
	dt = 0.01		# time step

	# --- Boundary Conditions (zero porosity for Neumann condition)
	u0, uL = 0.8, .1        # value at boundary, outside for Robin condition
	del0, delL = 1, 1       # porosity (reaction rate) for Robin condition

	# --- Variables
	x = np.linspace(0, L, Nx + 1)
	dx = L/Nx
	Nt = int(tf/dt)
	t = np.linspace(0, tf, Nt + 1)
	alpha = D/dx**2
	U = np.zeros( (Nx+1, Nt+1) )

	# --- Initial profile (Gaussian)
	uinit = np.exp(- (x - L/2)**2/(2*dx) )

	# --- Boundary conditions to pass to solver
	R = (2*dx/D) * np.array([del0, delL])
	B = (2*dx/D) * np.array([del0*u0, delL*uL])

	# --- Solve the IVP on the interior points
	soln = solve_ivp(de_rhs_NR, [0, tf], uinit, args=[alpha, R, B], dense_output=True)
	U = soln.sol(t)

	# --- Animation
	doMovie(x, t, U)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	diffusion_NR_via_MOL()
