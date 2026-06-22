#!/usr/bin/env python3
"""
Diffusion with Neumann or Robin BCs via Forward Euler

We approximate the solution of the Diffusion Equation (5.37) -- (5.38) with
using either Neumann or Robin boundary conditions using Forward Euler

Specifically, we find a numerical solution to the BVP

 u_t = D u_xx
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
# Forward Euler Solve
# =============================================================================
def FE_solve_NR(x, t, uinit, D, BC0, BCL):
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

    # --- Forward Euler Operator
    gam = D*dt/(dx**2)
    BFE = np.eye(Nx+1) + gam*D2

    # --- Nonhomogeneous BC
    zbc = np.zeros(Nx+1)
    zbc[0] = (2*dx/D) * (del0*U0)
    zbc[-1] = (2*dx/D) * (delL*UL)
    zbc = gam*zbc

    # --- Initialization
    U = np.zeros( (Nx+1, Nt+1) )
    U[:, 0] = uinit

    # --- Steps
    uk = uinit
    for kt in range(Nt):
        ukp1 = BFE@uk + zbc
        U[:, kt+1] = ukp1
        uk = ukp1
    return U
        
# =============================================================================
# Create Animation
# =============================================================================
def doMovie(x, t, U, ktskip):
    uinit = U[:,0]
    Nt = np.size(t) - 1

    # --- Initialization
    fig, ax = plt.subplots()
    p_update = ax.plot([], [], 'b', label='Time Evolution')[0]
    p_init = ax.plot(x, uinit, '--r', label='Initial Profile')
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

    ani = manimation.FuncAnimation(fig=fig, func=update, 
                frames=range(0, Nt+1, ktskip), interval=100)
    plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def FEuler_Diffusion_NR():
	# --- Global parameters
	L = 1
	D = .1

	# --- Boundary Conditions (zero porosity for Neumann condition)
	U0, UL = 1, 0		# value at boundary, outside for Robin condition
	del0, delL = 1, 1	# porosity (reaction rate) for Robin condition

	# --- Discretization
	Nx = 2**6		# number of spatial partitions
	Nt = 2**9		# number of temporal partitions
	dx = L/Nx		# spatial discretization
	dt = 0.001		# temporal discretization, dt < dx**2/(2D)
	tf = dt*Nt		# end time
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, tf, Nt+1) 

	# --- Initial profile, pass boundary conditions to solver
	u0_profile = np.exp( -(x - L/2)**2/(2*dx) )
	BC0 = np.array([del0, U0])
	BCL = np.array([delL, UL])

	# --- Solution and animation
	U = FE_solve_NR(x, t, u0_profile, D, BC0, BCL)
	doMovie(x, t, U, 2**3)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    FEuler_Diffusion_NR()
