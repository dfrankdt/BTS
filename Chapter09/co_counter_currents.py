#!/usr/bin/env python3
"""
Co and Counter Currents: We compute the transfer fraction as a function of
residence length

TO DO: Fix the legend it's trash
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def co_counter_currents():
	# --- Parameters
	L = np.linspace(0, 5, 2**9+1)
	rho_list = np.array([0.5, 2.0])
	
	# --- Initialize Plot
	fig, ax = plt.subplots()
	Co_plt = ax.plot([], [], '--k', label = 'Cocurrent')
	Cntr_plt = ax.plot([], [], '-k', label = 'Countercurrent')
	
	# --- Step through rho values
	for kp in range(len(rho_list)):
		rho = rho_list[kp]
		
		# --- cocurrent
		EL = np.exp(-L*(1 + 1/rho))
		c1L = (1 + rho*EL)/(1 + rho)
		c2L = rho*(1 - EL)/(1 + rho)
		
		# --- countercurrent
		EcL = np.exp(-L*(1 - 1/rho))
		c1cL = (1 - rho)*EcL/(EcL - rho)
		c2cL = rho*(EcL - 1)/(EcL - rho)
		
		# --- Do the plotting
		ax.plot(L, c2L, '--')
		ax.plot(L, c2cL, '-')
	
#	ax.legend(handles = [Co_plt, Cntr_plt], loc = 'upper right')
	ax.set(xlabel = r'Residence length, $dL/v_1$', ylabel = 'Transfer Friction')
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    co_counter_currents()
