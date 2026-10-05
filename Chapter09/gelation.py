#!/usr/bin/env python3
"""
Gelation: This script solves the pde

  dW/dt + d/dx (vW) = 0

where the initial profile is given by W0 and the velocity is prescribed
by a general form v = W/2 in order to give Burgers equation.

Produces: Animation illustrating the solution on the timeframe
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# Nonlinearity (Initial Profile)
# =============================================================================
def W(x, C0, k):
	Wofz = C0*k*(x - x**(k-1))
	return Wofz

# =============================================================================
# Velocity
# =============================================================================
def v(x, t, u):
	vel = 1/2*u
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
	ax.set(xlabel='z', ylabel='W(z, t)')
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
	return ani
	
# =============================================================================
# Create Snapshots
# =============================================================================
def doSnapShots(x, t, U):
	# --- Initialize data structures
	Nt = len(t) - 1
	
	fig, ax = plt.subplots()
	for kt in range(Nt+1):
		uk = U[:, kt]
		ax.plot(x, uk)
	ax.set(xlabel = 'z', ylabel = 'W(z, t)')
	return fig, ax

# =============================================================================
# Exact Solution
# =============================================================================
def doExact(z0, t, C0, k):
	# --- Initialize data structures
	Nt = len(t)-1
	Nx = len(x) - 1
	w = np.array( (Nx+1, Nt+1) )
	
	
	
# =============================================================================
# Main Simulation Function
# =============================================================================
def pde_upwind_MOL():
	# --- Discretizations
	L, Nx = 1, 100
	z = np.linspace(0, L, Nx+1)

	# --- Parameters
	C0 = 2/3
	k = 3

	# --- Initial Profile	
	w0 = W(z, C0, k)
	
	# --- Set the ODE
	tf, Nt = 1.5, 15
	soln = solve_ivp(de_rhs, [0, tf], w0, args=[z], dense_output=True)

	# --- Structure to produce visualization
	t = np.linspace(0, tf, Nt+1)
	W = soln.sol(t)
	ani = doMovie(z, t, W)
	fig, ax = doSnapShots(z, t, W)
	
	# --- Exact Solution
	W = doExact(
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    pde_upwind_MOL()
