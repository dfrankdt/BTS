#!/usr/bin/env python3
"""
RBC_plots: We create plots that look like figure 9.1 (more soon). 

NOTE probably add some titles
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Nonlinearity
# =============================================================================
def F(N):
	# --- Hill Function
	y = 1/(1 + N**7)
	return y

# =============================================================================
# Figure 9.1
# =============================================================================
def do_Fig_9_1(b_values, N):
	# --- Figure (a)
	fig1, ax1 = plt.subplots()
	ax1.plot(N, F(N))
	ax1.set(xlabel = 'U', ylabel = r'Production Rate $f(U)/A$')
	for kb in range(len(b_values)):
		b = b_values[kb]
		ax1.plot(N, b*N, label = rf'$\beta = ${b:1.2f}')
		ax1.legend(loc='upper right')
		ax1.set(xlim = (-0.1, 2.1), ylim = (-0.1, 1.1))
	
	# --- Figure (b) 
	fig2, ax2 = plt.subplots()
	ax2.plot(N/F(N), N)
	ax2.set(xlabel = 'XA', ylabel = r'$U_0$')
	ax2.set(xlim = (-0.1, 4.1), ylim = (-0.1, 1.5))
	
	return fig1, fig2

# =============================================================================
# Figure 9.2
# =============================================================================
def do_Fig_9_2():

# =============================================================================
# Main Simulation Function
# =============================================================================
def rbc_plots():
	# --- Global Parameters
	d = 70		# Time Delay
	X = 50		# RBC lifetime
	A = 1/50	# initial condition
	A = 0.09	
	
	# --- Steady state values
	b_values = np.array([0.2, 0.5, 0.8])
	N = np.linspace(0, 3, 2**8+1)

	# --- Do the plotting
	F91a, F91b = do_Fig_9_1(b_values, N)
	

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    rbc_plots()
