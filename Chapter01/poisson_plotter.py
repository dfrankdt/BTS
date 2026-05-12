#!/usr/bin/env python3
"""
Poisson Plotter: Create the Poisson distributions shown in Figure 1.7

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from math import factorial

# =============================================================================
# Main Simulation Function
# =============================================================================
def poisson_plotter():
	# --- Set parameters
	t = np.linspace(0, 10, 2**9+1)
	fig, ax = plt.subplots()
	ax.set(xlabel = r'$\alpha t$', ylabel = r'$p_j(t)$')

	for n_plot in range(5):
		j = n_plot
		f = t**j * np.exp(-t)/factorial(j)
		plot_label = [rf'$j = $ {j:1d}']
		ax.plot(t, f, label=plot_label)
	plt.legend()
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    poisson_plotter()
