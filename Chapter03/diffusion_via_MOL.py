#!/usr/bin/env python3
"""
Diffusion via Method of Lines

We simulate the diffusion equation with no flux boundary conditions using the
method of lines.  We use an odd number of cells with chemical concentration 1
in the middle cell and zero in the remaining cells.

The output of this code is an animation that illustrates diffusion of the 
chemical in the middle cell through the other cells.

This script is based on Diffusion_via_MOL.m
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# Right-hand Side of the Differential Equation
# =============================================================================
def de_rhs(t, u, alpha):
	N = len(u)-1
	
	# --- Second difference operator	
	D2u = np.zeros(N+1)
	D2u[1:N] = u[0:N-1] - 2*u[1:N] + u[2:N+1]
	
	# --- No flux boundary conditions
	D2u[0] = -2*u[0] + 2*u[1]
	D2u[-1] = 2*u[-2] - 2*u[-1]
	dudt = alpha*( D2u )
	return dudt

# =============================================================================
# Create Movie
# =============================================================================
def doMovie(x, t, U):
	# --- Initialize data structures
	Nt = len(t) - 1
	uinit = U[:,0]
	
	# --- Steady State
	uss = max(uinit)/(len(x) - 1)*np.ones(len(x))

	# --- Initialize movie
	fig, ax = plt.subplots()
	p_init = ax.plot(x, uinit, '.r', label='Initial Profile')
	p_ss = ax.plot(x, uss, '--y', label='Steady State')
	p_update = ax.plot([], [], '.b', label='Time Evolution')[0]
	ax.set(xlabel='n', ylabel='u(n, t)')
	ax.legend(loc='upper right')

	# --- Function to update the plot with the current frame
	def update(frame):
		tk = t[frame]
		u = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(u)
		ax.set(title=f'Time t = {tk:.2f} s')
		return(p_update)

	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt+1), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def diffusion_via_MOL():
	# --- Global parameters
	alpha = 0.1		# Time constant
	N = 41			# Spatial discretization (in this case, must be odd)
	tf = 100		# Final time for simulation
	Nt = 50			# Temporal discretization 

	# --- Numerical solution of the ODE
	uinit = np.zeros(N)
	uinit[int((N-1)/2)] = 1
	soln = solve_ivp(de_rhs, [0, tf], uinit, args=[alpha], dense_output=True)

	# --- Structure to produce visualization
	x = np.arange(N)
	t = np.linspace(0, tf, Nt+1)
	U = soln.sol(t)
	doMovie(x, t, U)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    diffusion_via_MOL()
