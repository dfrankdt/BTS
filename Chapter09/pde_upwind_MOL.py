#!/usr/bin/env python3
"""
Upwinding via Method of Lines

This is not yet correct
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
def de_rhs(t, u, x):
	Nx = len(x) - 1
	dx = x[-1]/Nx
	
	# --- Flux J = v*u, 
	xjmh = np.zeros(Nx+2)
	xjmh[:-1] = x - dx/2
	xjmh[-1] = x[-1] + dx/2
	
	ujmh = np.zeros(Nx+2)
	ujmh[1:-1] = (u[:-1] + u[1:])/2
	ujmh[0] = u[0]/2
	ujmh[-1] = u[-1]/2
	
	vjmh = xjmh
	ujm1 = np.zeros(Nx+2)
	ujm1[1:] = u
	
	uj = np.zeros(Nx+2)
	uj[:-1] = u
	Jmh = vjmh*ujm1*(vjmh>0) + vjmh*uj*(vjmh<0)	
	dJmh = (Jmh[:-1] - Jmh[1:])/dx
	du = dJmh
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
	p_exact = ax.plot([], [], 'g', label='Exact')[0]
	ax.set(xlabel='x', ylabel='u(x, t)')
	ax.legend(loc='upper right')

	# --- Function to update the plot with the current frame
	def update(frame):
		tk = t[frame]
		uk = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(uk)
		p_exact.set_xdata(x)
		p_exact.set_ydata( (x-tk)*(1-x+tk) )
		ax.set(title=f'Time t = {tk:.2f} s')
		return(p_update)

	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt+1), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def pde_upwind_MOL():
	# --- Discretizations
	Nx = 2**7
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
