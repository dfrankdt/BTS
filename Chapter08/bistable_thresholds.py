#!/usr/bin/env python3
"""
Thresholds in the Bistable Equation

We illustrate thresholds beyond which the Bistable Equation exhibits traveling
wave behavior.  We simulate the Bistable Equation for a range of length constants
and check whether the solution increases (i.e. becomes a traveling wave) or decays
from its initial profile.
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
#  Nonlinearity
# =============================================================================
def F(u, alpha):
	y = u * (1 - u) * (u - alpha)
	return y

# =============================================================================
# Crank-Nicolson Method
# =============================================================================
def doCN(x, t, uinit, D, alpha):
	""" 
	Crank-Nicolson to simulate the Bistable Equation

        u_t = Du_xx +  u (1 - u) (u-a)

	with no-flux boundary conditions.
	"""
	dx = x[1]-x[0]
	dt = t[1]-t[0]
	Nx = np.size(x) - 1
	Nt = np.size(t) - 1

	U = np.zeros( (Nx+1, Nt+1) )
	U[:,0] = uinit

	# -- Matrices for performing CN
	gam = D*dt/(2*dx**2)
	D2 = -2*np.eye(Nx+1)
	D2 = D2 + np.eye(Nx+1, k=1) + np.eye(Nx+1, k=-1)
	
	# -- No Flux BCs
	D2[0,0], D2[-1, -1] = -1, -1
	Acn = np.eye(Nx+1) - gam*D2
	Bcn = np.eye(Nx+1) + gam*D2
	
	# -- Initialization of CN method
	uk = uinit
	for kt in range(Nt):
		y = Bcn@uk + dt*F(uk, alpha)
		ukp1 = np.linalg.solve(Acn, y)
		U[:,kt+1] = ukp1
		uk = ukp1
	return U

# =============================================================================
# Traveling Wave Test
# =============================================================================
def doWaveCheck(x, t, D, alpha, a, lam):
	u0_profile = a*(1/np.cosh(x/lam))**2

	# --- Perform Crank-Nicolson, check
	U = doCN(x, t, u0_profile, D, alpha)
	z  = (U[0, -1] > alpha)
	return z

# =============================================================================
# Main Simulation Function
# =============================================================================
def bistable_thresholds():
	# --- Parameters
	L = 30
	D = 1
	tf = 50
	
	# --- Constants in initial profile to loop
	alpha_list = np.array([0.4,  0.25, 0.1])
	lam_list = np.linspace(0.1, 10, 50)
	
	# --- Discretization
	Nx = 2**6
	Nt = 2**8
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, tf, Nt+1)
	
	# --- Test through alpha values
	fig,ax = plt.subplots()
	ax.set(xlabel = r'length constant, $\lambda$', ylabel='Amplitude')
	ax.text(1, 0.1, 'subthreshold')#, size=12)
	ax.text(3.5, 0.8, 'superthreshold')#, size=12)
	ax.set(ylim=(0, 1))
	for kalpha in range(len(alpha_list)):
		critical_a = np.zeros(len(lam_list))
		alpha = alpha_list[kalpha]
		for klam in range(len(lam_list)):
			lam = lam_list[klam]
			mpb, mpa = 0, 1
			zoa = doWaveCheck(x, t, D, alpha, mpa, lam)

			# --- Bisection on the amplitude of the initial profile
			if (zoa == 1):
				for j in range(12):
					mpc = (mpa + mpb)/2
					zoc = doWaveCheck(x, t, D, alpha, mpc, lam)
					mpa = (1-zoc)*mpa + zoc*mpc
					mpb = zoc*mpb + (1-zoc)*mpc
				critical_a[klam] = mpc
			else:
				critical_a[klam] = 0
			# --- Remove superfluous values
			klam_threshold = np.nonzero(critical_a)
		ax.plot(lam_list[klam_threshold], critical_a[klam_threshold],
			label=rf'$\alpha$ = {alpha:1.2f}')
	ax.legend(loc='upper right')
	plt.show()
			
			# --- Check traveling wave, record
			

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	bistable_thresholds()
