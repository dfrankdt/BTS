#!/usr/bin/env python3
"""
PDE by Method of Characteristics: We solve

 du/dt + f(x, t, u) du/dx = g(x, t, u)
 
by first identifying a characteristic curve x = X0(t0), u(X0(t0), t0) = U0(t0) 
then using the chain rule to identify that along any such characteristic

 dX/dt = f(X, t, u)
 du/dt = g(X, t, u)

Note that the general conservation form

 du/dt + d/dx(vu) = 0
 
gives rise to

 du/dt + v du/dx = -(dv/dx) u
 
so that f(x, t, u) = v and g(x, t, u) = -(dv/dx) but the solution may be better
approximated by upwinding instead.


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
# Nonlinearities
# =============================================================================
"""
Note that
 - f = c, g = 0 gives rise to du/dt + c du/dx = 0
 - f = x, g = -1 gives rise to du/dt + d/dx (xu) = 0
 - f = 2u, g = 0 gives rise to du/dt + d/dx(u*u) = 0 (i.e. Burgers Equation)
all of which may be captured by upwinding.
"""
def f(x, t, u):
	#y = 1/2*np.ones(len(x))
	#y = x
	y = 2*u
	return y

def g(x, t, u):
	#y = np.zeros(len(x))
	#y = -u
	y = np.zeros(len(x))
	return y

# =============================================================================
# IVP RHS
# =============================================================================
def de_rhs(t, z):
	# --- Identify x, u from z
	Nx = int(len(z)/2 - 1)
	u = z[:Nx+1]
	x = z[Nx+1:]
	
	# --- Use nonlinearities to compute dx, du
	du = g(x, t, u)
	dx = f(x, t, u)
	
	# --- Construct dz
	dz = np.zeros(2*(Nx+1))
	dz[:Nx+1] = du
	dz[Nx+1:] = dx
	return dz

# =============================================================================
# Create Movie
# =============================================================================
def doMovie(X, t, U):
	# --- Initialize data structures
	Nt = len(t) - 1
	uinit = U[:,0]
	xinit = X[:,0]

	# --- Initialize movie
	fig, ax = plt.subplots()
	p_init = ax.plot(xinit, uinit, 'r', label='Initial Profile')
	p_update = ax.plot([], [], 'b', label='Time Evolution')[0]
	ax.set(xlabel='x', ylabel='u(x, t)')
	ax.legend(loc='upper right')

	# --- Function to update the plot with the current frame
	def update(frame):
		tk = t[frame]
		uk = U[:, frame]
		xk = X[:, frame]
		p_update.set_xdata(xk)
		p_update.set_ydata(uk)
		ax.set(title=f'Time t = {tk:.2f} s')
		return(p_update)

	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt+1), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def pde_by_MOC():
	# --- Discretizations
	L, Nx = 1, 2**5
	x0 = np.linspace(0, L, Nx+1)
	#u0 = x0*(1-x0)

	# --- Parameters and initial profile to solve Burgers equation in 9.5.2
	C0, k = 2/3, 3
	u0 = k*C0*(x0 - x0**(k-1))  
	
	# --- Set the ODEs
	tf, Nt = 1, 20
	z0 = np.zeros(2*(Nx+1))
	z0[:Nx+1] = u0
	z0[Nx+1:] = x0
	
	soln = solve_ivp(de_rhs, [0, tf], z0, dense_output = True)
	t = np.linspace(0, tf, Nt+1)
	z = soln.sol(t)
	U = z[:Nx+1,:]
	X = z[Nx+1:, :]
	doMovie(X, t, U)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    pde_by_MOC()
