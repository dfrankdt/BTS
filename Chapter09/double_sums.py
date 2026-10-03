#!/usr/bin/env python3
"""
Double Sums
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def double_sums():
	# --- Parameters
	N = 10
	
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 4.8))
	ax1.plot([1, N], [1, N], 'k')
	ax2.plot([1, N], [1, N], 'k')
	for n in range(1, N+1):
		ax1.plot(range(n,N+1), n*np.ones(N+1-n), '--.')
		ax2.plot(n*np.ones(n), range(1,n+1), '--.')
	ax1.set(xlim=(0.5, N+.5), ylim=(0.5, N+.5))
	ax1.set(xlabel = 'n', ylabel = 'i')
	ax2.set(xlim=(0.5, N+.5), ylim=(0.5, N+.5))
	ax2.set(xlabel = 'n', ylabel = 'i')
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    double_sums()
