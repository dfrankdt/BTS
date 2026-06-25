#!/usr/bin/env python3
"""
Hodgkin Huxley Upstroke

We illustrate bistability of the Hodgkin-Huxley equations by plotting the 
right-hand side.

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def HH_upstroke_rhs():
	# --- Parameters
	Vmin, Vmax = -20, 120
	gNabar = 120
	gKbar = 36
	gLbar = 0.3
	VNa = 115
	VK =  -12
	VL =  -10.5988
	h = 0.5961
	n = 0.3177

	# --- Functions
	def Ion_current(V):
		AM = 0.1*(25 - V)/(np.exp(0.1*(25 - V)) - 1)
		BM = 4*np.exp(-V/18)
		m = AM/(AM + BM)
		gNa = gNabar*m**3 * h
		gK = gKbar * n**4
		IL= gLbar*(V - VL)
		INa = gNa*(V - VNa)
		IK = gK*(V - VK)
		I = -IL -INa - IK
		return I
		
	# --- Variables
	N = 2**12
	V = np.linspace(Vmin, Vmax, N+1)
	Iion = Ion_current(V)

	# --- This identifies the zeros fairly well
	ndx = np.where(Iion[0:N]*Iion[1:N+1]<0)
	
	fig, ax = plt.subplots()
	ax.plot(V, Ion_current(V))
	ax.plot(V[ndx], Ion_current(V[ndx]), 'o')
	ax.set(xlabel='V (mV)', ylabel = r'-I$_{\text{ion}}$(V) ($\mu$A/cm$^2$)')

	print(f'Zeros of the Ion Current occur at V = {V[ndx][0]:1.2f}, {V[ndx][1]:1.2f}, and {V[ndx][2]:1.2f} mV')
	plt.show()


	
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    HH_upstroke_rhs()
