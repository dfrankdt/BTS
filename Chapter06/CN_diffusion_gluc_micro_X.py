#!/usr/bin/env python3
"""
Crank Nicolson scheme to simulate the diffusion with growth model where

	U + G -> 2U (rate alpha)
	
given in section 6.3.  Here we simulate the full non-scaled model on 0 < x < L
for 0 < t < tf.

Note: This script is based on CN_diffusion_gluc_micro_X.m
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# =============================================================================
#  Nonlinearity
# =============================================================================
def F(u, g, alpha):
	# -- Reaction term for the RD Equation
	y = alpha * u * g
	return y
# =============================================================================
# Crank-Nicolson Method
# =============================================================================
def doCN(x, t, u0_profile, g0_profile, Du, Dg, alpha):
	"""
	Crank-Nicolson to approximate the solution of the PDEs 
	
		u_t = Du u_xx + alpha u s
		s_t = Ds s_xx - alpha u s
			
	on the given spatio-temporal interval.
	
	Inputs:
		u0 (ndarray): Initial u profile 
		g0 (ndarray): Initial s profile
		x (ndarray): Spatial array
		t (ndarray): Temporal array
		Du, Dg (float): Diffusion coefficients
		alpha (float): Reaction rate
		
	Uses:
		dx, dt (float): spatial and temporal step sizes
		Nx, Nt (int): Number of spatial and temporal steps
		gam_u, gam_g (float): Coefficient for CN
		D2 (ndarray): Second difference matrix
		Acn_u, Bcn_u (ndarray): Left and right hand side for CN (u)
		Acn_g, Bcn_g (ndarray): Left and right hand side for CN (g)
		y_u, y_g (ndarray): Intermediate vectors for CN
		uk, ukp1 (ndarray): Current and next time steps (u)
		gk, gkp1 (ndarray): Current and next time steps (g)
		
	Outputs:
		U (ndarray): Approximate solution for u in space and time
		G (ndarray): Approximate solution for g in space and time
	"""
	dx = x[1]-x[0]
	dt = t[1]-t[0]
	Nx = np.size(x) - 1
	Nt = np.size(t) - 1

	# --- Second difference operator, adjusted for no flux BC
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	D2[0,0], D2[-1, -1] = -1, -1

	# -- CN Split
	gam_u = Du*dt/(dx**2)
	gam_g = Dg*dt/(dx**2)

	Acn_u = np.eye(Nx+1) - (gam_u/2)*D2
	Bcn_u = np.eye(Nx+1) + (gam_u/2)*D2
	Acn_g = np.eye(Nx+1) - (gam_g/2)*D2
	Bcn_g = np.eye(Nx+1) + (gam_g/2)*D2

	# -- Initialization of CN method
	U = np.zeros( (Nx+1, Nt+1) )
	U[:,0] = u0_profile
	
	G = np.zeros( (Nx+1, Nt+1) )
	G[:,0] = g0_profile

	# --- Steps
	uk, gk = u0_profile, g0_profile
	for kt in range(Nt):
		y_u = Bcn_u@uk + dt*F(uk, gk, alpha)
		y_g = Bcn_g@gk - dt*F(uk, gk, alpha)

		ukp1 = np.linalg.solve(Acn_u, y_u)
		gkp1 = np.linalg.solve(Acn_g, y_g)

		U[:,kt+1] = ukp1
		G[:,kt+1] = gkp1

		uk = ukp1
		gk = gkp1
	return U, G

# =============================================================================
# Create Movie
# =============================================================================
def doMovie(x, t, U, G):
	u0 = U[:,0]
	g0 = G[:,0]
	Nt = np.size(t) - 1

	# Initialize movie
	fig, (ax1, ax2) = plt.subplots(2, 1)
	p1_init = ax1.plot(x, u0, '--r', label='Initial Profile, u')
	p2_init = ax2.plot(x, g0, '--r', label='Initial Profile, g')
	p1_update = ax1.plot([], [], 'b', label='Time Evolution')[0]
	p2_update = ax2.plot([], [], 'g', label='Time Evolution')[0]
	ax1.set(ylabel=r'$u(x, t)$', ylim=(0,1))
	ax1.legend(loc='upper left')
	ax2.set(xlabel = 'x', ylabel = r'$g(x, t)$', ylim=(0,1))
	ax2.legend(loc='upper left')

	def update(frame):
	    tk = t[frame]
	    u = U[:, frame]
	    g = G[:, frame]
	    p1_update.set_xdata(x)
	    p1_update.set_ydata(u)
	    p2_update.set_xdata(x)
	    p2_update.set_ydata(g)
	    ax1.set(title=f'Time t = {tk:.2f} s')
	    return(p1_update, p2_update)
        
	ani = manimation.FuncAnimation(fig=fig, func=update, 
		frames=range(Nt+1), interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_diffusion_gluc_micro_X():

	# --- Global Parameters
	alpha = 0.1
	L = 200
	Du = 0.1
	Dg = 0.1

	# -- Spatial and Temporal Scales
	tend = 400
	Nt, Nx = 2**5, 2**7
	dt, dx = tend/Nt,  L/Nx
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, tend, Nt+1)

	# -- Initial profile for state variables
	u0_profile = 0.1*np.exp(- (x-L/2)**2/10)
	g0_profile = np.ones(Nx+1)/2 - u0_profile

	# -- Perform Crank-Nicolson
	U, G = doCN(x, t, u0_profile, g0_profile, Du, Dg, alpha)

	# -- Create Movie
	doMovie(x, t, U, G)	

# =============================================================================
# Execute the simulation if the script is run directly.
# =============================================================================
if __name__ == "__main__":
    CN_diffusion_gluc_micro_X()


