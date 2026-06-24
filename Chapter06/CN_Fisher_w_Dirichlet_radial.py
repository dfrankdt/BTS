#!/usr/bin/env python3
"""
CN Scheme to the 2D Fisher Equation (radial symmetry)

We solve the two-dimensional Fisher equation, leveraging radial symmetry, using
Crank-Nicolson

TO DO
 - Fix the second difference operator
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
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
def doCN(r, t, uinit, D, UR):
	"""
	Crank-Nicolson to simulate the Fisher equation with radial symmetry
	
	  u_t = D  nabla^2 u + ru (K - U)
	  
	with u(R, t) = 0 and du/dr(0, t) = 0
	"""
	dr = r[1] - r[0]
	dt = t[1] - t[0]
	Nr = len(r) - 1
	Nt = len(t) - 1
	
	# --- Interior discretizations
	rp = (r[2:Nr+1] + r[1:Nr])/2
	rm = (r[1:Nr] + r[0:Nr-1])/2
	
	# --- Second derivative operator
	D2p = rp[0:Nr-2]/r[1:Nr-1]
	D2m = rm[1:Nr-1]/r[2:Nr]
	D2 = -2*np.eye(Nr-1) + np.diag(D2p, k=1) + np.diag(D2m, k=-1)

	# --- Matrices for CN
	gam = D*dt/(dr**2)
	Acn = np.eye(Nr-1) - (gam/2)*D2
	Bcn = np.eye(Nr-1) + (gam/2)*D2
	
	# --- Initialization
	U = np.zeros( (Nr+1, Nt+1) )
	U[:, 0] = uinit
	
	# --- Steps
	uk = uinit
	for kt in range(Nt):
		# --- Incorporate boundary condition at r=0, r=R
		zbc = np.zeros(Nr-1)
		zbc[0] = gam * (rm[0]/r[1])*uk[0]
		zbc[-1] = gam * (rp[-1]/r[-2])*UR

		y = Bcn@uk[1:Nr] + zbc + dt*F(uk[1:Nr])
		ukp1 = np.zeros(Nr+1)
		ukp1[1:Nr] = np.linalg.solve(Acn, y)
		
		# --- Enforce boundary condition at r=0, r=R
		ukp1[0] = ukp1[1]
		ukp1[-1] = UR

		U[:, kt+1] = ukp1
		uk = ukp1
	return U
# =============================================================================
# Animation
# =============================================================================
def doMovie(r, t, U, ktskip):
	uinit = U[:,0]
	Nt = len(t) - 1

	# --- Structure for 2d Surface
	Nr = len(r) - 1
	theta = np.linspace(0, 2*np.pi, Nr+1)
	R, Theta = np.meshgrid(r, theta)
	X = R*np.cos(Theta)
	Y = R*np.sin(Theta)	
	
	def update(frame, zarray, plot):
		tk = t[frame]
		Uk = zarray[:, frame]
		plot[0].remove()
		plot[0] = ax.plot_surface(X, Y, Uk, cmap="copper")
		ax.set(title=f'Time t = {tk:.2f} s')
		return(plot)
		
	# --- Initialization
	fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
	plot = [ax.plot_surface(X, Y, uinit), cmap="copper", label='Surface']
	ax.set(xlabel='x', ylabel = 'y')
	ax.set(zlim=(0,1))
	
	ani = manimation.FuncAnimation(fig=fig, func=update
			frames=range(0, Nt+1), fargs=[U, plot], interval=100)
	plt.show()

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_Fisher_w_Dirichlet_radial():
	# --- Global parameters
	R = 5 # or R = 5
	D = 1
	UR = 0
	
	# --- Discretization
	Nr = 2**6
	Nt = 2**8
	tf = 8
	r = np.linspace(0, R, Nr+1)
	t = np.linspace(0, tf, Nt+1)
	
	# --- Initial profile
	u0_profile = 0.2*np.exp(-r**2)
	U = doCN(r, t, u0_profile, D, UR)
	doMovie(r, t, U, 2**3)


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_Fisher_w_Dirichlet_radial()
