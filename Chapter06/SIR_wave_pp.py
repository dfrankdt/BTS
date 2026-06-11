#!/usr/bin/env python3
"""
SIR Phase Plane Wave

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
# Equilibrium Analysis
# =============================================================================
def get_s0(eta):
	"""
	The nonzero equilibrium satisfies 1 + eta log(s) - s = 0 provided eta < 1.
	"""
	sn = 0.1
	check = 1
	while check > 1e-8:
		f = 1 + eta*np.log(sn) - sn
		fp = eta/sn - 1
		snp1 = sn - f/fp
		check = np.abs(snp1 - sn)
		sn = snp1
	return sn

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(t, y, c, eta):
	s, u = y
	ds = s*u
	du = c*(1 + eta*np.log(s) - s - u)
	dy = np.array([ds, du])
	return dy

# =============================================================================
# Main Simulation Function
# =============================================================================
def SIR_wave_pp():
	# --- Parameters
	eta = 0.5
	s0 = get_s0(eta)
	c_list = np.array([1/0.8, 1/1.8])
	
	
	# --- Trajectories
	tend = 100
	t = np.linspace(0, tend, 2**8+1)
	yinit = [s0, 0.001]
	fig, ax = plt.subplots()
	for c in c_list:
		soln = solve_ivp(de_rhs, [0, tend], yinit, args=[c, eta], dense_output=True)
		S, U = soln.sol(t)
		ax.plot(S, U, label=f'1/c = {1/c:1.1f}')

	# --- Plotting
	snc = np.linspace(0.1, 1.1, 2**8+1)
	unc = 1 + eta*np.log(snc) - snc
	ax.plot([s0, 1], [0, 0], 'ok')
	ax.plot(snc, unc, '--y')
	ax.plot([0, 1.2], [0, 0], '--g')
	ax.set(xlim=(0, 1.2), ylim=(-0.05, 0.2))
	ax.legend()


	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    SIR_wave_pp()
