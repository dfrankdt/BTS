#!/usr/bin/env python3
"""
CN Scheme to solve Fisher's equation

We approximate the solution to Fisher's equation on 0 < x < L
(dimensionless) with zero Dirichlet boundary data using two
different values of L.

Produces: As in Figure 6.7
 - Animation showing time evolution of traveling wave (or not)
 - Figure displaying time snapshots of traveling wave (or not)
 
This script is based on CN_Fisher_w_Dirichlet.m 
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# =============================================================================
#  Nonlinearity
# =============================================================================
def F(u):
	y = u * (1-u)
	return y

# =============================================================================
#  CN Scheme
# =============================================================================
def doCN(x, t, uinit, D, BC):
	dx = x[1]-x[0]
	dt = t[1]-t[0]
	Nx = len(x) - 1
	Nt = len(t) - 1
	
	# --- Second difference operator
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)

	# --- CN split, adjusted for Dirichlet BC
	gam = D*dt/(dx**2)
	Acn = np.eye(Nx+1) - (gam/2)*D2
	Bcn = np.eye(Nx+1) + (gam/2)*D2
	
	Acn[0,:] = np.zeros(Nx+1)
	Acn[0, 0] = 1
	Acn[-1,:] = np.zeros(Nx+1)
	Acn[-1, -1] = 1

	# --- Initialization of CN method
	U = np.zeros( (Nx+1, Nt+1) )
	U[:,0] = uinit

	# --- Steps
	uk = uinit
	for kt in range(Nt):
		y = Bcn@uk + dt*F(uk)
		y[0], y[-1] = BC
		ukp1 = np.linalg.solve(Acn, y)
		U[:,kt+1] = ukp1
		uk = ukp1
	return U
	
# =============================================================================
#  Create Animation
# =============================================================================
def doMovie(x, t, U, ktskip):
	Nt = len(t) - 1
	u0 = U[:,0]
        
	# --- Initialize animation
	fig, ax = plt.subplots()
	p_update = ax.plot([], [], 'b', label='Time Evolution')[0]
	p_init = ax.plot(x, u0, '--r', label='Initial Profile')
	ax.set(ylim=(0,1))
	ax.set(xlabel='x', ylabel='u(x, t)')
	ax.legend(loc='upper left')

	# --- Update frame
	def update(frame):
		tk = t[frame]
		u = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(u)
		ax.set(title=f'Time t = {tk:.2f} s')
		ax.legend(loc='upper left')
		return(p_update)
        
	ani = manimation.FuncAnimation(fig=fig, func=update, 
			frames=range(0, Nt+1, ktskip), interval=100)

	plt.show()

# =============================================================================
#  Create Snapshots
# =============================================================================
def doSnapshots(x, t, U):
	# --- Initialize data structures, including snapshot frequency
	ksnap = 2**4
	dt = t[ksnap]
	Nt = len(t) - 1
	u0 = U[:,0]

	# --- Initialize animation
	fig, ax = plt.subplots()
	ax.plot(x, u0, '--r', label='Initial Profile')
	ax.set(ylim=(0,1))
	ax.set(xlabel=r'$\xi$', ylabel=r'v($\xi$, t)')
	
	for kt in range(ksnap, Nt+1, ksnap):
		u = U[:, kt]
		ax.plot(x, u)
	ax.set(title = f'Snapshots every t = {dt:1.2f} s')
	plt.show()
# =============================================================================
#  Main Simulation Function
# =============================================================================
def CN_Fisher_w_Dirichlet():
	# --- Global parameters: Choose a domain length to obtain a traveling wave
	D = 1	# Diffusion coefficient (WLOG)
##	L = 2	# Length 1: Below threshold for traveling wave
	L = 5	# Length 2: Above threshold for traveling wave
	
	# --- Boundary Conditions
	U0, UL = 0, 0

	# --- Discretization
	Nx = 2**6
	Nt = 2**9
	dx = L/Nx
	dt = 0.01
	tf = dt*Nt
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, tf, Nt+1)

	# --- Choose an initial profile, pass boundary conditions to solver
	BC = np.array([U0, UL])

	# --- Profile 1: Falls to traveling wave
	u0_profile = np.tanh(x/(L/25)) * np.tanh(-(x - L)/(L/25))
	# --- Profile 2: Rises to traveling wave
##	u0_profile = (1 + np.tanh( (x - L/2)/(L/25) ))*(1 + np.tanh(-(x - L/2)/(L/25)))/5
	
	# --- Solution, animation, and snapshots
	U = doCN(x, t, u0_profile, D, BC)
	doMovie(x, t, U, 2**2)
	doSnapshots(x, t, U)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	CN_Fisher_w_Dirichlet()
