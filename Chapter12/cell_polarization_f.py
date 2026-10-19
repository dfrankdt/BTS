#!/usr/bin/env python3
"""
Cell Polarization: We find a numerical solution to the reaction-diffusion equation
with a slowly varying stimulus S(x, t).

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation
rng = np.random.default_rng()

# =============================================================================
# Parameters
# =============================================================================
def getPars():
	Du = 0.001	# Diffusion for u
	Dv = 1		# Diffusion for v
	delta = 1	# Removal rate (WLOG)
	K = 1
	gamma = 0.5
	return Du, Dv, delta, K, gamma	

# =============================================================================
# Stimulus
# =============================================================================
def stim(x, t):
	# --- Parameters local to the stimulation: see eqn (12.43)
	ampmx = 0.1
	amprate = 0.005
	ampzero = np.pi/amprate
	
	# --- Amplitude
	signal = ( (1 - np.cos(2*np.pi*x))/2 )**10
	kamp = (t < ampzero)*ampmx*np.sin(amprate*t)**2
	z = kamp*signal
	return z

# =============================================================================
# Nonlinearity
# =============================================================================
def f(u, v, S):
	Du, Dv, delta, K, gamma = getPars()
	z = S + gamma*u**2/(K + u**2)
	return z*v - delta*u

# =============================================================================
# RD Solver
# =============================================================================
def doSolve(dt, Nt, Nt_inter, x, uinit, vinit):
	Nx = len(x) - 1
	dx = x[1] - x[0]
	Du, Dv, delta, K, gamma = getPars()
	
	# --- Initialize solution
	t = np.zeros( (Nt+1) )
	U = np.zeros( (Nx+1, Nt+1) )
	U[:, 0] = uinit
	V = np.zeros( (Nx+1, Nt+1) )
	V[:, 0] = vinit
	
	# --- CN Structure
	gamu = Du*dt/(2*dx**2)
	gamv = Dv*dt/(2*dx**2)
	D2 = -2*np.eye(Nx + 1)
	D2 = D2 + np.eye(Nx + 1, k=1) + np.eye(Nx + 1, k = -1)
	# --- No Flux
	D2[0, 0], D2[-1, -1] = -1, -1
	
	Acnu = np.eye(Nx+1) - gamu*D2
	Acnv = np.eye(Nx+1) - gamv*D2
	Bcnu = np.eye(Nx+1) + gamu*D2
	Bcnv = np.eye(Nx+1) + gamv*D2
	
	# --- Double sum to include intermediate steps
	uk, vk = uinit, vinit
	tk = 0
	for kt in range(Nt):
		for kt_inter in range(Nt_inter):
			S = stim(x, tk)
			yu = Bcnu@uk + dt*f(uk, vk, S)
			yv = Bcnv@vk - dt*f(uk, vk, S)
			ukp1 = np.linalg.solve(Acnu, yu)
			vkp1 = np.linalg.solve(Acnv, yv)
			
			uk = ukp1
			vk = vkp1
			tk += dt
		U[:, kt+1] = ukp1
		V[:, kt+1] = vkp1
		t[kt+1] = tk
	return t, U, V

# =============================================================================
# Animation
# =============================================================================
def doMovie(x, t, U):
	# --- Initialize data structures
	Nt = len(t) - 1
	uinit = U[:,0]
	Umax = max(U[int((len(x) - 1)/2), :])

	# --- Initialize movie
	fig, ax = plt.subplots()
	p_init = ax.plot(x, uinit, '--r', label='Initial Profile')
	p_update = ax.plot([], [], 'b', label='Time Evolution')[0]
	ax.set(xlabel='x', ylabel='u(x, t)', ylim = (0, Umax))
	ax.legend(loc='upper right')

	def update(frame):
		tk = t[frame]
		u = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(u)
		ax.set(title=f'Time t = {tk:.2f} s')
		return(p_update)
	ani = manimation.FuncAnimation(fig=fig, func=update,
			frames=range(0, Nt+1), interval=100)

	return ani

# =============================================================================
# Main Simulation Function
# =============================================================================
def cell_polarization_f():
	# --- Discretization
	L = 1		# Length (WLOG)
	Nx = 2**7
	Nt = 2**8
	Nt_inter = 2**5
	dt = 0.1
	x = np.linspace(0, L, Nx+1)
	
	# --- Initial profile
	Vt = 4.2
	u0 = np.zeros(len(x))
	v0 = Vt*np.ones(len(x))
	t, U, V = doSolve(dt, Nt, Nt_inter, x, u0, v0)
	
	ani = doMovie(x, t, U)
	
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    cell_polarization_f()
