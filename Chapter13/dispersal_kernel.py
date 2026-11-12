#!/usr/bin/env python3
"""
Invasion: We construct the dispersal kernel given by a symmetric, exponentially
decreasing function

Produces
 - Figure 1: Graph of the dispersal kernel as a function of scaled space

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def dispersal_kernel():
	# --- Variables
	L, Nx = 5, 2**8
	x = np.linspace(-L, L, 2*Nx+1)
	q = 1/2 * np.exp(-np.abs(x))
	
	# --- Dispersal Kernel
	fig, ax = plt.subplots()
	ax.plot(x, q)
	ax.set(xlabel = r'$\alpha$x', ylabel = r'K(x)/$\alpha$')
	
	plt.show()
	


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    dispersal_kernel()
