#!/usr/bin/env python3
"""
Invasion: We construct the profile of invasion given by a piecewise defined
function

Produces:
 - Figure 1: Graph of the invasion profile as a function of scaled space

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def invasion_profile():
	# --- Variables
	L, Nx = 5, 2**8
	x = np.linspace(-L, L, 2*Nx+1)
	fm = (1 - 1/2*np.exp(x))
	fp = 1/2*np.exp(-x)
	f = fm*( x<0 ) + fp * (x >= 0)
	
	# --- Dispersal Kernel
	fig, ax = plt.subplots()
	ax.plot(x, f)
	ax.set(xlabel = r'$\alpha$x', ylabel = r'$u_{n+1}$/$A_0$')
	
	plt.show()
	


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    invasion_profile()
