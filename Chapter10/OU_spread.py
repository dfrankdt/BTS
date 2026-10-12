#!/usr/bin/env python3
"""
OU Spread: We illustrate the spread of an OU Epidemic via SIR.

Produces:
 - Figure 1: Graph of sqrt(x)exp(-x) as a function of x, see Figure 10.4
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def OU_spread():
	# --- Variables
	Nx = 2**8
	x = np.linspace(0, 5, Nx+1)
	f = np.sqrt(x)*np.exp(-x)
	
	# --- Plotting
	fig, ax = plt.subplots()
	ax.plot(x, f)
	ax.set(xlabel = r'$\Lambda$', ylabel = r'$f(\Lambda)$')
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    OU_spread()
