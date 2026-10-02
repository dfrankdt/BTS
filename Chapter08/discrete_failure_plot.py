#!/usr/bin/env python3
"""
Discrete Failure: We plot the bifurcation diagram separating propagation from
propagation failure for the discrete bistable equation, described by eqn (8.65).

Produces
 - Figure 1: Bifurcation diagram, see Figure 8.11
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def discrete_failure_plot():
	# --- Variables
	a = np.linspace(0, 1, 2**8)
	d = a*(1-a)/(2*a-1)**2
	
	fig, ax = plt.subplots()
	ax.plot(a, d)
	ax.set(xlim = (0, 1), ylim = (0, 10))
	ax.set(xlabel = r'$\alpha$', ylabel = 'd')
	ax.text(0.2, 5, 'Propagation', ha = 'center')
	ax.text(0.8, 5, 'Propagation', ha = 'center')
	ax.text(0.5, 1, 'Propagation Failure', ha = 'center')
	
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    discrete_failure_plot()
