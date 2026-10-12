#!/usr/bin/env python3
"""
OU Process: We simulate an OU process by finding a numerical approximation to
the solution of the reaction-diffusion equation given by (10.17).

Produces:
 - Figure 1: Snapshots of the solution of the OU process as a function of time
 - Figure 2: 
 - Figure 3:
 
TO DO: Update figures 2 and 3
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Numerical Solve
# =============================================================================
def doSolve(x, uinit, D, kappa):
	# --- Spatial Discretization
	dx = x[1] - x[0]
	Nx = len(x) - 1
	
	# --- Spatial midpoints
	xm = np.zeros(Nx+2)
	xm[1:Nx+1] = (x[0:Nx] + x[1:Nx+1])/2
	xm[0] = x[0] - dx
	xm[-1] = x[-1] + dx
	
	# --- Temporal Discretization
	Nt = 2**5
	Nt_steps = 2**6
	dt = 0.9*dx/max(kappa*xm)
	t = np.zeros(Nt+1)
	
	# --- State variable
	U = np.zeros( (Nx+1, Nt+1) )
	U[:, 0] = uinit

	# --- Matrices for CN
	gam = D*dt/(2*dx**2)
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	# --- No Flux BCs
	D2[0, 0], D2[-1, -1] = -1, -1
	
	Acn = np.eye(Nx+1) - gam*D2
	Bcn = np.eye(Nx+1) + gam*D2
	
	# --- Double sum to include intermediate steps
	uk = uinit
	tk = 0
	for kt in range(Nt):
		for kt_inter in range(Nt_steps):
			# --- Flux for upwinding
			v = -kappa*(xm)
			kx_up = (v>0)
			upos = np.zeros(Nx+2)
			upos[1:] = uk
			uminus = np.zeros(Nx+2)
			uminus[:-1] = uk
			
			F = v * ( (v>0)*upos + (1 - (v>0))*uminus )
			dF = (F[1:] - F[:-1])/dx
			y = Bcn@uk - dt*dF
			ukp1 = np.linalg.solve(Acn, y)
			uk = ukp1
			tk += dt
		U[:, kt+1] = ukp1
		t[kt+1] = tk
	return t, U
			
# =============================================================================
# Main Simulation Function
# =============================================================================
def OU_process():
	# --- Parameters
	d = 1
	K = 1
	
	# --- Discretizations
	L = 4
	Nx = 2**8
	x = np.linspace(-L, L, Nx+1)
	
	# --- Initialization
	x0 = -2
	u0 = np.exp(-(x-x0)**2/0.1)
	
	# --- Numerical solution
	t, U = doSolve(x, u0, d, K)
	
	# --- Plotting
	fig, ax = plt.subplots()
	ax.plot(x, u0, '--')
	ax.plot(x, U[:, 1:])
	ax.set(xlabel = r'$\xi$', ylabel =r'$p(\xi, \tau)$')
	
	plt.show()
	
	# --- Structure to check the total integral
	int_u0 = np.trapezoid(u0, x)
	int_U = np.trapezoid(U, x, axis=0)/int_u0
	mj = np.argmax(U, 0)
	mxU = np.max(U, 0)

	fig2, ax2 = plt.subplots()
	ax2.plot(t, x[mj])
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    OU_process()
