#!/usr/bin/env python3
"""
Fire-Diffuse-Fire Model:

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
# Main Simulation Function
# =============================================================================
def fdf_plots():
	# --- Parameters
	n = 1000	# large
	dt = np.linspace(20/2**8, 20, 2**8)
	blist = np.array([0.25, 0.05, 0.01])
	colorlist = np.array(['b', 'g', 'y'])
	
	fig, ax = plt.subplots()
	for kb in range(len(blist)):
		b = blist[kb]
		c = colorlist[kb]
		theta = 0
		for kn in range(1, n+1):
			dtheta = 1/np.sqrt(4*np.pi*kn*dt)*np.exp(-kn/(4*dt) - b**2*kn*dt)
			theta += dtheta
		knmax = np.argmax(theta)
		ax.plot(dt[:knmax], theta[:knmax], color=[c],
				dt[knmax:], theta[knmax:], )
	plt.show()
			
	

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    fdf_plots()
