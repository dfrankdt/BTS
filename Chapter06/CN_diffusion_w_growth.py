#!/usr/bin/env python3
"""
Crank-Nicholson Reaction/Diffusion

We approximate the solution of the reaction diffusion equation

 u_t = D u_xx + sigma alpha u

using a Crank-Nicolson scheme.  As in the text, we incorporate the
parameter sigma to determine whether the equation exhibits growth
or decay

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
def F(u, sigma, alpha):
        y = sigma * alpha *u
        return y

# =============================================================================
# CN Scheme
# =============================================================================
def doCN(u0, x, t, D, sigma, alpha):
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
		y = Bcn@uk + dt*F(uk, sigma, alpha)
		ukp1 = np.linalg.solve(Acn, y)
		uk = ukp1
		U[:,kt+1] = ukp1
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
# Create Snapshots (Figure 6.2)
# =============================================================================
def doSnapshots(x, t, Ug, Ud):
        Nt = len(t) - 1
        kts = range(int(Nt/8), Nt+1, int(Nt/8))
        
        fig1, ax1 = plt.subplots()
        ax1.plot(x, Ug[:,0], '--b')
        ax1.plot(x, Ug[:, kts])
        ax1.set(xlabel = 'x', ylabel = 'u(x, t)', title='Growth')
        ax1.set(ylim=(0, 1))
        

        fig2, ax2 = plt.subplots()
        ax2.plot(x, Ud[:, 0], '--b')
        ax2.plot(x, Ud[:, kts])
        ax2.set(xlabel = 'x', ylabel = 'u(x, t)', title='Decay')
        ax2.set(ylim=(0, 1))
        plt.show()

# =============================================================================
# Create Log Plot (Figure 6.3)
# =============================================================================
def doLogPlot(x, t, Ug, Ud, U0):
        Nx = len(x) - 1

        UgRatio = Ug/U0
        UdRatio = Ud/U0

        fig3, ax3 = plt.subplots()
        ax3.semilogy(t, UgRatio[int(Nx/2),:], 
                         label = r'Growth, $\alpha \sigma = 0.2$')
        ax3.semilogy(t, UdRatio[int(Nx/2),:],
                         label = r'Growth, $\alpha \sigma = -0.2$')
        ax3.set(xlabel = 't', ylabel = r'$\log(U(20, t)/U_0(20, t))$')
        ax3.legend(loc='upper left')
        plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_diffusion_w_growth():
        # --- Global parameters
        L = 40
        D = 1
        alpha = 0.2

        Nt, Nx = 2**6, 2**8
        dt, dx = 0.2,  L/Nx

        # --- Set spatial and temporal scales, initial profile
        x = np.linspace(0, L, Nx+1)
        t = np.linspace(0, dt*Nt, Nt+1)
        u0 = 0.2*np.exp(- (x-L/2)**2/(L/4))

        # --- Do the solving
        print(f'We approximate the solution on 0 < x < {x[-1]:2.1f} for 0 < t < {t[-1]:2.1f}')
        
        
        sigma = -1
        U_decay = doCN(u0, x, t, D, sigma, alpha)
        sigma = 1
        U_growth = doCN(u0, x, t, D, sigma, alpha)
        sigma = 0
        U_diffusion = doCN(u0, x, t, D, sigma, alpha)

        # --- Animation or Snapshots, depending
        doMovie(x, t, U_decay)
        doMovie(x, t, U_growth)
        doSnapshots(x, t, U_growth, U_decay)
        doLogPlot(x, t, U_growth, U_decay, U_diffusion)




# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
        CN_diffusion_w_growth()



