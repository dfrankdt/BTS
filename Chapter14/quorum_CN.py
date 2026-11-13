#!/usr/bin/env python3
"""
CN Analysis of Quorum Sensing: 

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# Parameters
# =============================================================================
def doPars():
	#--- Parameters
	s0 = 0.05		# Basal production rate
	K = 1			# Production rate Kd
	g = 1.5			# Scale factor
	d = 1			# Transmembrane permeability
	alph = 0.0375   # Extracellular decay
	return s0, K, g, d, alph

# =============================================================================
# Nonlinearity
# =============================================================================
def f(u):
	s0, K, g, d, alph = doPars()
	y = s0 + g * u**2/(K**2 + u**2)
	return y

# =============================================================================
# CN Solve
# =============================================================================
def doSolve(x, u0, v0, rho, D):
	# --- Parameters
	s0, K, g, d, alph = doPars()

	# --- Discretization
	Nx = len(x) - 1
	dx = x[1] - x[0]
	dt = dx
	Nt = 2**4
	Nt_inter = 2**10
	t = np.zeros(Nt+1)
	
	# --- State Variables
	U = np.zeros( (Nx+1, Nt+1) )
	V = np.zeros( (Nx+1, Nt+1) )
	U[:, 0] = u0
	V[:, 0] = v0
	
	# --- Matrices for CN
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	# --- No Flux
	D2[0, 0], D2[-1, -1] = -1, -1
	
	Acn = np.eye(Nx+1) - (D*dt/(2*dx**2))*D2
	Bcn = np.eye(Nx+1) + (D*dt/(2*dx**2))*D2
	
	# --- Initialization of CN
	uk = u0
	vk = v0
	tk = 0
	for kt in range(Nt):
		for kt_inter in range(Nt_inter):
			yu = uk + dt*(f(uk) + d*(vk - uk))
			yv = Bcn@vk + dt*(-alph*vk + rho/(1-rho)*d*(uk - vk))
			
			ukp1 = yu
			vkp1 = np.linalg.solve(Acn, yv)
			
			uk = ukp1
			vk = vkp1
			tk += dt
		# --- Record
		U[:, kt+1] = uk
		V[:, kt+1] = vk
		t[kt+1] = tk
		
	return t, U, V
	
# =============================================================================
# Animation
# =============================================================================
def doMovie(x, t, U, V):
	# --- Initialize
	Nt = len(t) - 1
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (6.4, 9.6))
	u_init = ax1.plot(x, U[:,0], '--r')
	u_update = ax1.plot([], [], 'b')[0]	
	v_update = ax2.plot([], [], 'b')[0]
	ax1.set(xlabel = 'x', ylabel = 'Intracellular Autoinducer, u(x, t)', 
		xlim = (0, 10), ylim=(0, 2))
	ax2.set(xlabel = 'x', ylabel = 'Intracellular Autoinducer, v(x, t)', 
		xlim = (0, 10), ylim=(0, 2))

	# --- Animation
	def update(frame):
		uk = U[:, frame]
		vk = V[:, frame]
		u_update.set_xdata(x)
		u_update.set_ydata(uk)
		v_update.set_xdata(x)
		v_update.set_ydata(vk)
		return(u_update, v_update)
		
	ani = manimation.FuncAnimation(fig = fig, func = update, frames = range(1, Nt+1),
		interval = 100)
	return ani

# =============================================================================
# Main Simulation Function
# =============================================================================
def quorum_CN():
	# --- Discretization
	D = 0.001
	L, Nx = 10, 2**7
	x = np.linspace(0, L, Nx+1)
	
	# --- Initialization and solve
	rho = 0.02
	u0 = 2*(1/np.cosh(2*(x - L/2))**2)
	v0 = np.zeros(Nx+1)
	t, U, V = doSolve(x, u0, v0, rho, D)
	
	# ---
	fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (6.4, 9.6))
	ax1.plot(x, U[:, 1:])
	ax1.set(xlabel = 'x', ylabel = 'Intracellular Autoinducer, u(x, t)', ylim=(0, 2))
	ax2.plot(x, V[:, 1:])
	ax2.set(xlabel = 'x', ylabel = 'Intracellular Autoinducer, v(x, t)', ylim=(0, 2))
	

	ani = doMovie(x, t, U, V)	
	plt.show()
	
	


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	quorum_CN()
