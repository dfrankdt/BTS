#!/usr/bin/env python3
"""
Diffusion with Neumann or Robin BCs via Forward Euler

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
#    zbc = gam*zbc

    # --- Initialization
    U = np.zeros( (Nx+1, Nt+1) )
    U[:, 0] = uinit

    # --- Steps
    uk = uinit
    for kt in range(Nt):
        ukp1 = BFE@uk + dt*zbc
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
    del0, U0 = 0, 0
    delL, UL = 0, 0

    # --- Discretization
    Nt, Nx = 2**8, 2**6
    dx = L/Nx
    dt = 0.0001      # need dt < dx**2/(4D)
    tf = dt*Nt
    x = np.linspace(0, L, Nx+1)
    t = np.linspace(0, tf, Nt+1) 

    # --- Boundary Conditions, initial profile
    BC0 = np.array([del0, U0])
    BCL = np.array([delL, UL])
    u0_profile = np.exp( -(x - L/2)**2/(2*dx) )

    # --- Solution and animation
    U = FE_solve_NR(x, t, u0_profile, D, BC0, BCL)
    doMovie(x, t, U, 2**3)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    FEuler_Diffusion_NR()
