#!/usr/bin/env python3
"""
Fisher Equation with Neumann or Robin conditions via Method of Lines

We approximate the solution of the Fisher equation using either Neumann or Robin boundary conditions with the method of lines

Specifically, we find a numerical solution to the BVP

 u_t = D u_xx + r u(K - u)
 Du_x = delta_0(u - U0) at x=0
 -Du_x = deltaL(u - UL) at x=L

where the nonhomogenous boundary conditions are given as in equation (5.23). In this case, the capital letters indicate given values at the boundary, but outside the domain.
"""


# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
#  Nonlinearity
# =============================================================================
def F(u, r, K):
	y = r*u * (K - u)
	return y
# =============================================================================
# DE RHS -- Neumann or Robin Condition
# =============================================================================
def de_rhs_NR(t, u, dx, D, r, K, del0, U0, delL, UL):
	N = len(u)-1
	alpha = D/dx**2

	D2u = np.zeros(N+1)
	D2u[1:N] = u[0:N-1] - 2*u[1:N] + u[2:N+1]
	D2u[0] = (-2 - 2*dx/D*del0)*u[0] + 2*u[1] + 2*dx/D*del0*U0 
	D2u[-1] = 2*u[-2] + (-2 - 2*dx/D*delL)*u[-1] + +2*dx/D*delL*UL
	dudt = alpha*D2u + F(u, r, K)
	return dudt

# =============================================================================
# Animation
# =============================================================================
def doMovie(x, t, U, ktskip):
	# --- Initialize data structures
	Nt = len(t) - 1
	uinit = U[:,0]
	
	# --- Initialize movie
	fig, ax = plt.subplots()
	p_init = ax.plot(x, uinit, '-r', label='Initial Profile')
	p_update = ax.plot([], [], '-b', label='Time Evolution')[0]
	ax.set(xlabel='x', ylabel='u(x, t)')
	ax.set(ylim=(0,1))
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
def MOL_Fisher_NR():
	# --- Global parameters
	L = 100			# length of domain
	tf = 40
	D = 1			# diffusion coefficient
	K = 1			# Carrying Capacity
	r = 1			# base growth rate

	# --- Boundary Conditions (zero porosity for Neumann condition)
	U0, UL = 1, 0		# value at boundary, outside for Robin condition
	del0, delL = 0, 0	# porosity (reaction rate) for Robin condition
	
	# --- Discretization
	Nx = 2**10		# number of spatial partitions
	Nt = 2**7		# number of temporal partitions
	dx = L/Nx		# spatial discretization
	dt = tf/Nt		# temporal discretization
#	tf = dt*Nt		# end time
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, tf, Nt+1)

 	# --- Initial profile, pass boundary conditions to solver
	IVP_arg = [dx, D, r, K, del0, U0, delL, UL]
	u0_profile = 0.002*np.exp(- (x - L/2)**2/(2*dx) )

	# --- Solve the IVP
	soln = solve_ivp(de_rhs_NR, [0, tf], u0_profile, args=IVP_arg, dense_output=True)
	U = soln.sol(t)

	# --- Animation
	doMovie(x, t, U, 2**0)


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    MOL_Fisher_NR()
