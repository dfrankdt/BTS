#!/usr/bin/env python3
"""
Diffusion with Dirichlet BCs via Crank-Nicolson

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# Crank-Nicolson Solve
# =============================================================================
def CN_solve(x, t, uinit, D, BC):
	dx = x[1] - x[0]
	dt = t[1] - t[0]
	Nx = len(x) - 1
	Nt = len(t) - 1

	# --- Second difference operator
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	
	# --- CN split, adjusted for BC
	gam = D*dt/(dx**2)
	Acn = np.eye(Nx+1) - gam/2*D2
	Bcn = np.eye(Nx+1) + gam/2*D2
	Acn[0,] = np.zeros(Nx+1)
	Acn[0, 0] = 1
	Acn[-1,] = np.zeros(Nx+1)
	Acn[-1, -1] = 1
	
	# --- Initialization
	U = np.zeros( (Nx+1, Nt+1) )
	U[:, 0] = uinit
	
	uk = uinit
	for kt in range(Nt):
		y = Bcn@uk
		y[0], y[-1] = BC
		ukp1 = np.linalg.solve(Acn, y)
		U[:, kt+1] = ukp1
		uk = ukp1
	return U
# =============================================================================
# Create Animation
# =============================================================================
def doMovie(x, t, U):
	u0 = U[:,0]
	Nt = np.size(t) - 1

	# --- Initialization
	fig, ax = plt.subplots()
	p_update = ax.plot([], [], 'b', label='Time Evolution')[0]
	p_init = ax.plot(x, u0, '--r', label='Initial Profile')
	ax.set(ylim=(0,1))
	ax.set(xlabel='x', ylabel='u(x, t)')
	ax.legend(loc='upper left')

	# --- Update
	def update(frame):
		tk = t[frame]
		u = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(u)
		ax.set(title=f'Time t = {tk:.2f} s')
		ax.legend(loc='upper left')
		return(p_update)
        
	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt+1), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_diffusion_D():
	# --- Global parameters
	L = 1
	D = .1
	u0, uL = 0, 0

	# --- Discretization
	Nt, Nx = 250, 2**6
	dx = L/Nx
	dt = 0.01#dx**2/(4*D)
	tf = dt*Nt
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, tf, Nt+1) 
	
	# --- Boundary conditions, initial profile
	BC = np.array([u0, uL])
	u0_profile = np.exp( -(x - L/2)**2/(2*dx) )

	# --- Solution and animation
	U = CN_solve(x, t, u0_profile, D, BC)
	doMovie(x, t, U)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_diffusion_D()
