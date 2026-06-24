#!/usr/bin/env python3
"""
Bessel Function

Plot the Bessel function of the first kind, J0(x). We access the SciPy special
functional library.

Figures Produced: J0(r) as a function of 0 <= r <= 25.

This script is basd on J0plot.m 

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import j0

# =============================================================================
# Main Simulation Function
# =============================================================================
def J0plot():
	# --- Global parameters
	Nx, R = 2**8, 25
	r = np.linspace(0, R, Nx+1)

	# --- Do the plotting	
	fig, ax = plt.subplots()
	ax.plot(r, j0(r))
	ax.plot([0, R], [0, 0], '--r')
	ax.set(ylim=(-1, 1), xlabel = r'$r$', ylabel=r'$J_0(r)$')
	
	plt.show()
	


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    J0plot()
