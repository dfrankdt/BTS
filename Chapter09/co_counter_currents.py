#!/usr/bin/env python3
"""
Co and Counter Currents: We compute the transfer fraction as a function of
residence length

Produces
 - Figure 1: Illustration of transfer fraction as a function of residence length
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
	colorlist = ['b', 'g']
	
	# --- Initialize Plot
	fig, ax = plt.subplots()
	Cntr_plt, = ax.plot([], [], '-k', label = 'Countercurrent')
	Co_plt, = ax.plot([], [], '--k', label = 'Cocurrent')
	ax.legend(handles = [Cntr_plt, Co_plt], loc = 'upper left')
	ax.set(xlabel = r'Residence length, $dL/v_1$', ylabel = 'Transfer Fraction')
	
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
		ax.plot(L, c2L, '--', color = colorlist[kp])
		ax.plot(L, c2cL, '-', color = colorlist[kp])
	
	# --- Annnotate	
	ax.text(2.5, 0.375, rf'$\rho = ${rho_list[0]:1.2f}')
	ax.text(2.5, 0.725, rf'$\rho = ${rho_list[1]:1.2f}')
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    co_counter_currents()
