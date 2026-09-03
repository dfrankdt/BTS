#!/usr/bin/env python3
"""
CN Scheme to the 2D Fisher Equation (radial symmetry)

We solve the two-dimensional Fisher equation, leveraging radial symmetry, using
Crank-Nicolson

TO DO: Does the initial profile really matter? Probably choose one.
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
def F(u):
	y = u*(1 - u)
	return y

# =============================================================================
# CN Solve
# =============================================================================
def doCN(r, t, uinit, D):
	"""
	Crank-Nicolson to simulate the Fisher equation with radial symmetry
	
	  u_t = D  nabla^2 u + ru (K - U)
	  
	with u(R, t) = 0 and du/dr(0, t) = 0
	"""
	dr = r[1] - r[0]
	dt = t[1] - t[0]
	Nr = len(r) - 1
	Nt = len(t) - 1
	
	# --- Second derivative operator
	cm = np.linspace(1, Nr-1, Nr-1)
	cp = np.linspace(0, Nr-2, Nr-1)
	cp[0] = 1/2		# Adjust for zero flux at r = 0
	D2p = np.diag(1 + 1/(2*cp), k=1)
	D2m = np.diag(1 - 1/(2*cm), k=-1)

	D2 = D2m - 2*np.eye(Nr) + D2p
	
	# --- Matrices for CN
	gam = D*dt/(dr**2)
	Acn = np.eye(Nr) - (gam/2)*D2
	Bcn = np.eye(Nr) + (gam/2)*D2
	
	# --- Initialization
	U = np.zeros( (Nr+1, Nt+1) )
	U[:, 0] = uinit
	
	# --- Steps
	uk = uinit[0:Nr]
	for kt in range(Nt):
		# --- Solve the linear system, CN on Laplacian, forward Euler on F(u)
		y = Bcn@uk + dt*F(uk)
		ukp1 = np.linalg.solve(Acn, y)
		# --- Record solution and step forward
		U[0:Nr, kt+1] = ukp1
		uk = ukp1
	return U
# =============================================================================
# Animations
# =============================================================================
def do3dMovie(r, t, U, ktskip):
	uinit = U[:,0]
	Nt = len(t) - 1

	# --- Structure for 2d Surface
	Nr = len(r) - 1
	theta = np.linspace(0, 2*np.pi, Nr+1)
	R, Theta = np.meshgrid(r, theta)
	X = R*np.cos(Theta)
	Y = R*np.sin(Theta)
	Uinit, Ones = np.meshgrid(uinit, np.ones(Nr+1))
	
	# --- Initialization
	fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
	plot = [ax.plot_surface(X, Y, Uinit, cmap="copper", label='Surface')]
	ax.set(xlabel='x', ylabel = 'y')
	ax.set(zlim=(0,1))

	# --- Animation update
	def update(frame, zarray, plot):
		tk = t[frame]
		Uk, Ones = np.meshgrid(zarray[:, frame], np.ones(Nr+1))
		plot[0].remove()
		plot[0] = ax.plot_surface(X, Y, Uk, cmap="copper")
		ax.set(title=f'Time t = {tk:.2f} s')
		return(plot)
		
	
	ani = manimation.FuncAnimation(fig=fig, func=update,
			frames=range(0, Nt+1, ktskip), fargs=[U, plot], interval=100)
	return ani

def do2dMovie(r, t, U, ktskip):
	uinit = U[:,0]
	Nt = len(t) - 1

	# --- Structure for curve
	Nr = len(r) - 1
	
	# --- Initialization
	fig, ax = plt.subplots()
	p_init = ax.plot(r, uinit, '--r', label = 'Initial Profile')
	p_update = ax.plot([], [], 'b', label = 'Time Evolution')[0]
	ax.set(xlabel='r', ylabel = 'u(r, t)')
	ax.set(ylim = (0, 1))
	ax.legend(loc = 'upper right')

	# --- Animation update
	def update(frame):
		tk = t[frame]
		uk = U[:, frame]
		p_update.set_xdata(r)
		p_update.set_ydata(uk)
		ax.set(title=f'Time t = {tk:.2f} s')
		return(p_update)
		
	
	ani = manimation.FuncAnimation(fig=fig, func=update,
			frames=range(0, Nt+1, ktskip), interval=100)
	return ani
	
# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_Fisher_w_Dirichlet_radial():
	# --- Global parameters
	R = 4 # or R = 10
	D = 1
	UR = 0
	
	# --- Discretization
	Nr = 2**6
	Nt = 2**8
	tf = 15
	r = np.linspace(0, R, Nr+1)
	t = np.linspace(0, tf, Nt+1)
	
	# --- Initial profile
	u0_profile = 0.2*np.exp(-(r-R/2)**2)
	u0_profile = 1 - np.tanh( (r - 0.8)/(1/20))
	u0_profile = 0.01 * (1 - np.tanh((r**2)/(1/20)))
	u0_profile = 0.2*np.exp(-r**2)
	u0_profile[-1] = UR
	U = doCN(r, t, u0_profile, D)
	ani2d = do2dMovie(r, t, U, 2**3)
	ani3d = do3dMovie(r, t, U, 2**3)
	
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_Fisher_w_Dirichlet_radial()
