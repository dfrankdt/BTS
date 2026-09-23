#!/usr/bin/env python3
"""
Fire-Diffuse-Fire Model:

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def fdf_plots():
	# --- Parameters: For this simulation a logarithmically spacing is appropriate
	n = 1000
	dt = np.logspace(np.log(25/2**8)/np.log(10), np.log(25)/np.log(10), 2**8)
	blist = np.array([0.01, 0.05, 0.25])

	# --- Initialize Plot
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize= (12.8, 4.8) )
	colorlist = np.array(['b', 'g', 'r'])

	# --- Step through beta values	
	for kb in range(len(blist)):
		b = blist[kb]
		theta = 0
		for kn in range(1, n+1):
			dtheta = 1/np.sqrt(4*np.pi*kn*dt)*np.exp(-kn/(4*dt) - b**2*kn*dt)
			theta += dtheta
		# --- Capture maximum value: we are interested in the increasing curve
		knmax = np.argmax(theta)
		
		# --- Theta as a function of dt
		ax1.plot(dt[:knmax], theta[:knmax], '-', color = colorlist[kb],
			label = rf'$\beta = ${b:1.2f}')
		ax2.plot(theta[:knmax], 1/dt[:knmax], '-', color = colorlist[kb], 
			label = rf'$\beta = ${b:1.2f}')
			
		# --- 1/dt as a function of Theta
		ax1.plot(dt[knmax:], theta[knmax:], '--', color = colorlist[kb])
		ax2.plot(theta[knmax:], 1/dt[knmax:], '--', color = colorlist[kb])

	# --- Labeling
	ax1.set(xlabel = r'$\delta_\tau$', ylabel = r'$\theta = g_\beta(\delta_\tau)$')
	ax2.set(xlabel = r'$\theta = c^*L/\sigma$', ylabel = r'$1/\delta_\tau$')
	ax1.legend(loc = 'upper right')
	ax2.legend(loc = 'upper right')
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    fdf_plots()
