#!/usr/bin/env python3
"""
Diffusion with Neumann or Robin BCs via Crank-Nicolson

We find a numerical solution to the BVP

 u_t = D u_xx
 Du_x = delta_0(u - U0) at x=0
 -Du_x = deltaL(u - UL) at x=L

where the nonhomogenous boundary conditions are given as in equation (5.23). In this case, the capital letters indicate given values at the boundary, but outside the domain.

TO DO: Adjust figsize and zlim
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.animation as manimation

# =============================================================================
# Crank-Nicolson
# =============================================================================
def doCN(uinit, x, t, D, BC0, BCL):
	dx = x[1] - x[0]
	dt = t[1] - t[0]
	Nx = len(x) - 1
	Nt = len(t) - 1
	
	del0, U0 = BC0
	delL, UL = BCL
	
	# --- Second difference operator, adjusted for BC
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	D2[0, 0] = -2*dx/D*del0 - 2
	D2[0, 1] = 2
	D2[-1, -2] = 2
	D2[-1, -1] = -2*dx/D*delL - 2
	
	# --- Nonhomogenous BC
	zbc = np.zeros(Nx+1)
	zbc[0] = (2*dx/D) * (del0*U0)
	zbc[-1] = (2*dx/D) * (delL*UL)
	
	# --- CN Split
	gam = D*dt/(dx**2)
	Acn = np.eye(Nx+1) - gam/2*D2
	Bcn = np.eye(Nx+1) + gam/2*D2
	
	# --- Initialization
	U = np.zeros( (Nx+1, Nt+1) )
	U[:,0] = uinit
	
	uk = uinit
	for kt in range(Nt):
		y = Bcn@uk + dt*zbc
		ukp1 = np.linalg.solve(Acn, y)
		U[:,kt+1] = ukp1
		uk = ukp1
	return U

# =============================================================================
# Animation
# =============================================================================
def doMovie(x, t, U):
	u0 = U[:,0]
	Nt = len(t) - 1
	
	# --- Initialization
	fig, ax = plt.subplots()
	p_init = ax.plot(x, u0, '--r', label = 'Initial Profile')
	p_update = ax.plot([], [], 'b', label = 'Time Evolution')[0]
	ax.set(xlabel = 'x', ylabel = 'u(x, t)')
	ax.legend(loc = 'upper left')
	
	# --- Update
	def update(frame):
		tk = t[frame]
		u = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(u)
		ax.set(title=f'Time t = {tk:.3f} s')
		return(p_update)

	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(0, Nt+1), interval=100)

	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_Diffusion_NR():
	# --- Parameters
	D = 1
	L = 1
	tend = 0.04
		
	# --- Discretizations
	Nt, Nx = 250, 2**8
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, tend, Nt+1)
	dt, dx = tend/Nt, L/Nx

	# --- Boundary Conditions
	del0, u0 = 0, 0	# For Robin -Du_x + del0 u = del0 u0 at x = 0
	delL, uL = 0, 0 # For Robin  Du_x + delL u = delL uL at x = L
	BC0 = np.array([del0, u0])
	BCL = np.array([delL, uL])
	
	# --- Initialization and surface u(x, t)
	uinit = 1.5*np.cos(np.pi*(x - L/2))
	U_plt_a = doCN(uinit, x, t, D, BC0, BCL)
	
	uinit = 1.5*np.cos(np.pi*(x - L/2)) + 1.5*np.cos(20*np.pi*x)
	U_plt_b = doCN(uinit, x, t, D, BC0, BCL)
	
	T, X = np.meshgrid(t, x)
	fig1, (axa, axb) = plt.subplots(1, 2, 
		figsize=(12.8, 4.8),
		subplot_kw={"projection": "3d"})
	axa.plot_surface(X, T, U_plt_a, cmap="coolwarm")
	axa.set(zlim=(-1, 2))
	axb.plot_surface(X, T, U_plt_b, cmap="coolwarm")
	axb.set(zlim=(-1, 2))
	plt.show()


	# --- Initialization and animation
	uinit = 1.5*np.exp(-(x - L/2)**2/(2*dx))
	U = doCN(uinit, x, t, D, BC0, BCL)
	doMovie(x, t, U)

	
	
	
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_Diffusion_NR()
