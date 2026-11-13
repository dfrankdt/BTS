#!/usr/bin/env python3
"""
Quorum Sensing Plots: We illustrate the existence of hysteresis in steady 
solutions of a model describing quorum sensing. 

Produces
 - Figure 1: Autoinducer u as a function of p = u/f(u), showing that gamma
   sufficiently large produces possibility of a biphasic response

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def quorum_plots():
	# --- Parameters
	K = 1
	s0 = 0.1
	gam_values = np.array([0.5, 1, 2.5])
	
	# --- Graphing
	u = np.linspace(0, 10, 2**8+1)
	fig, ax = plt.subplots()
	for kg in range(len(gam_values)):
		g = gam_values[kg]
		f = s0 + g*u**2/(K**2 + u**2)
		p = u/f
		ax.plot(p, u, label = rf'$\gamma$ = {g:1.1f}')
	ax.legend(loc = 'upper left')
	ax.set(xlabel = 'p', ylabel = 'u', xlim=(0, 5), ylim=(0, 5))
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    quorum_plots()
