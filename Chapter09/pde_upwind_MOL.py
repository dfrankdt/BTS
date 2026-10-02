#!/usr/bin/env python3
"""
Upwinding via Method of Lines: This script solves the pde

  du/dt + d/dx (vu) = 0

where the initial profile is given by u0 and the velocity is prescribed
by a general form v = v(x, t, u). Below we can prescribe the initial profile
and the velocity.


"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# Velocity
# =============================================================================
def v(x, t, u):
	vel = x
	return vel

# =============================================================================
# Right-hand Side of the Differential Equation
# =============================================================================
def de_rhs(t, u, x):
	# --- u, x have Nx + 1 entries. 
	Nx = len(x) - 1
	dx = x[1] - x[0]
	
	# --- Need the flux J = v*u, but first need half step values for x and u 
	xjmh = np.zeros(Nx+2)
	xjmh[:-1] = x - dx/2
	xjmh[-1] = x[-1] + dx/2

	um = np.zeros(Nx+2)
	um[:-1] = u
	up = np.zeros(Nx+2)
	up[1:] = u
	
	ujmh = (um+up)/2

	# --- Call the velocity and get the flux, be sure to upwind.	
	vjmh = v(xjmh, t, ujmh)
	Jmh = vjmh*( (vjmh>0)*up + (vjmh<0)*um )
	
	# --- Now get the difference
	dJ = Jmh[1:] - Jmh[:-1]
	du = -dJ/dx
	
	return du
	
# =============================================================================
# Create Movie
# =============================================================================
def doMovie(x, t, U):
	# --- Initialize data structures
	Nt = len(t) - 1
	uinit = U[:,0]

	# --- Initialize movie
	fig, ax = plt.subplots()
	p_init = ax.plot(x, uinit, 'r', label='Initial Profile')
	p_update = ax.plot([], [], 'b', label='Time Evolution')[0]
	ax.set(xlabel='x', ylabel='u(x, t)')
	ax.legend(loc='upper right')

	# --- Function to update the plot with the current frame
	def update(frame):
		tk = t[frame]
		uk = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(uk)
		ax.set(title=f'Time t = {tk:.2f} s')
		return(p_update)

	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt+1), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def pde_upwind_MOL():
	# --- Discretizations
	Nx = 2**5
	L = 1
	x = np.linspace(0, L, Nx+1)
	u0 = x*(1-x)
	
	# --- Set the ODE
	Nt = 20
	tf = 1
	soln = solve_ivp(de_rhs, [0, tf], u0, args=[x], dense_output=True)

	# --- Structure to produce visualization
	t = np.linspace(0, tf, Nt+1)
	U = soln.sol(t)
	doMovie(x, t, U)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    pde_upwind_MOL()
