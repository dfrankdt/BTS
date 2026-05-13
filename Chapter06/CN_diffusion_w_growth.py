#!/usr/bin/env python3
"""
Crank-Nicholson Reaction/Diffusion

We approximate the solution of the reaction diffusion equation

 u_t = D u_xx + alpha u

using a Crank-Nicolson scheme. The differential equation is linear,
so the growth term may be incorporated into the differential operator.

This script is based on CN_diffusion_w_growth.m
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# =============================================================================
# Reaction Term
# =============================================================================
def F(u):
        alpha = .2
        k = 1
        y = alpha*k*u
        return y

# =============================================================================
# CN Scheme
# =============================================================================
def doCN(u0, x, t, D):
	dx = x[1]-x[0]
	dt = t[1]-t[0]
	Nx = np.size(x) - 1
	Nt = np.size(t) - 1

	gam = D*dt/(2*dx**2)
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)

	Acn = np.eye(Nx+1) - gam*D2
	Bcn = np.eye(Nx+1) + gam*D2

	U = np.zeros( (Nx+1, Nt+1) )
	U[:,0] = u0

	uk = u0
	for kt in range(Nt):
		y = Bcn@uk + dt*F(uk)
		ukp1 = np.linalg.solve(Acn, y)
		uk = ukp1
		U[:,kt+1] = ukp1
	return U

# =============================================================================
# Create Animation
# =============================================================================
def doMovie(x, t, U):
        Nt = np.size(t) - 1
        # --- Initialization
        u0 = U[:,0]
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
def CN_diffusion_w_growth():
        # --- Global parameters
        L = 40
        D = 1

        Nt, Nx = 2**5, 2**8
        dt, dx = 0.2,  L/Nx

        x = np.linspace(0, L, Nx+1)
        t = np.linspace(0, dt*Nt, Nt+1)

        u0 = 0.2*np.exp(- (x-L/2)**2/(L/4))
        U = doCN(u0, x, t, D)
        doMovie(x, t, U)


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
        CN_diffusion_w_growth()



