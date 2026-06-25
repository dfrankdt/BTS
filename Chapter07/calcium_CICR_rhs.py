#!/usr/bin/env python3
"""
CICR: Calcium induced calcium release

Bistability in CICR

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def calcium_CICR_rhs():
	# --- Global Parameters
	Vp = 0.9
	kf = 1.1
	ct = 2
	K1 = 0.1
	K2 = 0.08
	gm = 5.5
	Kp = 0.1
	
	# --- Flux functions
	def Jserca(c):
		z = Vp*c**2/(Kp**2 + c**2)
		return z
	
	def JIPR(c):
		P = p**3/(p + K1)**3
		P0 = P*c**3/(c + K2)**3
		ce = gm*(ct - c)
		z = kf*P0*(ce - c)
		return z

	# --- Varying parameters	
	p_list = np.array([0.12, 0.15, 0.2])
	c = np.linspace(0, ct, 2**6+1)
	
	fig, ax = plt.subplots()
	ax.plot(c, np.zeros(len(c)), '--k')
	for kp in range(len(p_list)):
		p = p_list[kp]
		Jtot = 	JIPR(c) - Jserca(c)
		ax.plot(c, Jtot, label=f'p = {p:.2f}')
	ax.set(xlabel=r'Ca$^{++}$ ($\mu$M)', ylabel = r'$J_\text{total}$')
	ax.set(ylim=(-0.8, 0.8))
	ax.legend()
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    calcium_CICR_rhs()
