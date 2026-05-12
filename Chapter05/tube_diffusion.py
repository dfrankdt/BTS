#!/usr/bin/env python3
"""
Tube Diffusion

We create the plots in Figure 5.1 of the text. We access the SciPy error function
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

# =============================================================================
# Main Simulation Function
# =============================================================================
def tube_diffusion():
	# --- Global parameters
	Nz, zf = 2**8, 4
	u0 = 1

	# --- z and u(z) for plotting
	z = np.linspace(zf/Nz, zf, Nz)		# We avoid z = 0 for Figure 2
	u = u0 * (1 - erf(z))
	
	# --- Do the plotting
	fig1, ax1 = plt.subplots()
	ax1.plot(z, u)
	ax1.set(xlabel = 'z', ylabel = 'u(z)')
	
	fig2, ax2 = plt.subplots()
	ax2.plot(np.sqrt(1/z), u)
	ax2.set(xlabel = r'$z^{-1/2}$', ylabel = 'u(z)')
	
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    tube_diffusion()
