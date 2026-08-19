#!/usr/bin/env python3
"""
RBC_plots: We create plots that look like figure 9.1 (more soon)

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def rbc_plots():
	# --- Parameters
	b_values = np.array([0.2, 0.5, 0.8])
	N = np.linspace(0, 3, 2**8+1)
	F = 1/(1 + N**7)

	# --- Figure (a)
	fig, ax1 = plt.subplots()
	ax1.plot(N, F)
	ax1.set(xlabel = 'U', ylabel = r'Production Rate $f(U)/A$')
	for kb in range(len(b_values)):
		b = b_values[kb]
		ax1.plot(N, b*N, label = rf'$\beta = ${b:1.2f}')
		ax1.legend()
		ax1.set(ylim = (-0.1, 1.1))
	
	# --- Figure (b) 
	fig, ax2 = plt.subplots()
	ax2.plot(N/F, N)
	ax2.set(xlabel = 'XA', ylabel = r'$U_0$')
	ax2.set(xlim = (-0.1, 4.1), ylim = (-0.1, 1.5))	
	plt.show()


	do_prod_rate()
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    rbc_plots()
