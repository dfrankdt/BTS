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
# DE RHS
# =============================================================================
def de_rhs(t, y, c, eta):
	u, s = y
	ds = 1/c*s*u
	du = c*(1 + eta*np.log(s) - s - u)
	dy = np.array([ds, du])
	return dy

# =============================================================================
# Main Simulation Function
# =============================================================================
def SIR_wave_pp():
	# --- Parameters
	eta = 0.47
	
	# --- Plotting
	snc = np.linspace(0.1, 1.1, 2**8+1)
	unc = 1 + eta*np.log(snc) - snc
	fig, ax = plt.subplots()
	ax.plot([0, 1], [0, 0], 'ok')
	ax.plot(snc, unc, '--r')	
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    SIR_wave_pp()
