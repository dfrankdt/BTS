#!/usr/bin/env python3
"""
Crank-Nicholson Scheme to the 2d Fisher Equation (Dirichlet BCs)

We solve the two-dimensional Fisher equation (6.74) with homogeneous Dirichlet 
conditions. As noted, there is one parameter

u_t = u_xx + u_yy + mu u(1-u)

We note that this particular code examines a square grid with N partitions in
the x and y directions. Given the homogeneous Dirichlet data, we solve only at 
the (N-1) x (N-1) interior points.

TO DO:
 - Possibly incorporate the initial profile in the movie
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# =============================================================================
# Nonlinearity
# =============================================================================
def F(u, mu):
	y = mu*u * (1-u)
	return y

# =============================================================================
# CN Scheme
# =============================================================================
def doCN(x, y, t, uinit, mu):
	N = len(x) - 1
	Nt = len(t) - 1
	h = x[1] - x[0]
	dt = t[1] - t[0]
	
	# --- Kronecker Product for Laplacian
	D2 = -2*np.eye(N-1) + np.eye(N-1, k=1) + np.eye(N-1, k=-1)
	I = np.eye(N-1)
	L = np.kron(D2, I) + np.kron(I, D2)
	
	# --- CN Splitting
	gam = dt/h**2
	ACN = np.eye( (N-1)**2 ) - (gam/2)*L
	BCN = np.eye( (N-1)**2 ) + (gam/2)*L
	
	# --- Initialization
	Uinterior = np.zeros( (N-1, N-1, Nt+1) )
	Uinterior[:,:,0] = uinit
	
	# --- Steps
	uk = np.reshape(uinit, (N-1)**2)
	for kt in range(Nt):
		y = BCN@uk
		ukp1 = np.linalg.solve(ACN, y) + dt*F(uk, mu)
		Uinterior[:, :, kt+1] = np.reshape(ukp1, (N-1, N-1))
		uk = ukp1
		
	# --- Embed solution in zero boundary conditions
	U = np.zeros( (N+1, N+1, Nt+1) )
	U[1:N, 1:N, :] = Uinterior
	return U

# =============================================================================
# Create Movie
# =============================================================================
def doMovie(x, y, t, U):
	# --- Initialize data structures
	Nt = np.size(t) - 1
	uinit = U[:,:,0]
	X, Y = np.meshgrid(x, y)
	

	def update(frame, zarray, plot):
	    tk = t[frame]
	    Uk = zarray[:, :, frame]
	    plot[0].remove()
	    plot[0] = ax.plot_surface(X, Y, Uk, cmap="copper")
	    ax.set(title=f'Time t = {tk:.2f} s')
	    return(plot)
        
	# --- Initialize movie
	fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
	plot = [ax.plot_surface(x, y, uinit, cmap="winter", label='Surface')]
	ax.set(xlabel='x', ylabel='y')
	ax.set(zlim=(0, 1))

	ani = manimation.FuncAnimation(fig=fig, func=update, 
			frames=range(0, Nt+1), fargs=[U, plot], interval=100)
	plt.show()


# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_Fisher_w_Dirichlet_2D():
	# --- Global parameters
	mu = 22 # not exactly sure about threshold, yes traveling wave mu = 22
	
	# --- Discretizations
	N = 2**5
	x = np.linspace(0, 1, N+1)
	y = np.linspace(0, 1, N+1)
	Nt = 2**6
	tf = 3.5
	t = np.linspace(0, tf, Nt+1)
	
	# --- Initial profile, surface Gaussian
	[X0, Y0] = np.meshgrid(x[1:N], y[1:N])
	u0_profile = np.exp(- (X0-0.5)**2/(0.01))*np.exp(-(Y0-0.5)**2/(0.01))/10
	
	fig0, ax0 = plt.subplots(subplot_kw={"projection": "3d"})
	ax0.plot_surface(X0, Y0, u0_profile, cmap="winter")
	plt.show()
	U = doCN(x, y, t, u0_profile, mu)
	doMovie(x, y, t, U)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_Fisher_w_Dirichlet_2D()
