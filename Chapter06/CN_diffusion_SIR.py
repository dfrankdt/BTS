#!/usr/bin/env python3
"""
CN Scheme to the RD Equation with SIR Dynamics

We solve the system

 s_t = - su 
 u_t = su + eta u + u_xx
 
described in equations (6.100) - (6.101) where s represents the dimensionless
susceptible population and u represents the dimensionless infected population.
The relevant dimensionless parameter is eta

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
def F(s, u, alpha):
	# --- Reaction term for the RD Equation
	y = alpha * s * u
	return y

# =============================================================================
# Crank-Nicolson Method
# =============================================================================
def doCN(x, t, sinit, uinit, D, eta, r):
	dx = x[1] - x[0]
	dt = t[1] - t[0]
	Nx = len(x) - 1
	Nt = len(t) - 1
	
	# --- Second difference operator, adjusted for no flux
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	D2[0,0] = -1
	D2[-1, -1] = -1
	
	# --- CN Split
	gam = D*dt/(dx**2)
	Acn = np.eye(Nx+1) - (gam/2)*D2
	Bcn = np.eye(Nx+1) + (gam/2)*D2

	# --- Initialization	
	S = np.zeros( (Nx+1, Nt+1) )
	U = np.zeros( (Nx+1, Nt+1) )
	S[:,0] = sinit
	U[:,0] = uinit

	# --- Steps	
	sk = sinit
	uk = uinit
	for kt in range(Nt):
		y = Bcn@uk + dt * (F(sk, uk, r) - eta*uk)
		skp1 = sk - dt*(F(sk, uk, r))
		ukp1 = np.linalg.solve(Acn, y)
		S[:, kt+1] = skp1
		U[:, kt+1] = ukp1
		sk = skp1
		uk = ukp1
	return S, U
		

# =============================================================================
# Create Movie
# =============================================================================
def doMovie(x, t, S, U, ktskip):
	s0 = S[:,0]
	u0 = U[:,0]
	Nt = len(t) - 1
	
	# --- Initialization
	fig, (ax1, ax2) = plt.subplots(2, 1)
	p1_init = ax1.plot(x, u0, '--b', label='Initial Profile')
	p2_init = ax2.plot(x, s0, '--b', label='Initial Profile')
	p1_update = ax1.plot([], [], 'r', label='Time Evolution')[0]
	p2_update = ax2.plot([], [], 'g', label='Time Evolution')[0]
	ax1.set(ylabel=r'$u(\xi, \tau)$', ylim=(0,1))
	ax1.legend(loc='upper left')
	ax2.set(xlabel = r'$\xi$', ylabel = r'$\sigma(\xi, \tau)$', ylim=(0,1))
	ax2.legend(loc='upper left')

	# --- Update
	def update(frame):
		tk = t[frame]
		u = U[:, frame]
		s = S[:, frame]
		p1_update.set_xdata(x)
		p1_update.set_ydata(u)
		p2_update.set_xdata(x)
		p2_update.set_ydata(s)
		ax1.set(title=f'Time t = {tk:.2f} s')
		return(p1_update, p2_update)
        
	ani = manimation.FuncAnimation(fig=fig, func=update, 
			frames=range(0,Nt+1, ktskip), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_diffusion_SIR():
	# --- Parameters
	D = 1
	r = 1
	eta = 0.5
	L = 80

	# --- Spatial and Temporal Scales
	tend = 60
	Nt, Nx = 2**8, 2**7
	dt, dx = tend/Nt,  L/Nx
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, tend, Nt+1)

	# --- Initial Profile
	u0_profile = 0.1 / np.cosh(2*(x - L/2))**2
	s0_profile = 1 - u0_profile
	
	# --- Solve and animate
	S, U = doCN(x, t, s0_profile, u0_profile, D, eta, r)
	doMovie(x, t, S, U, 2**3)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_diffusion_SIR()
