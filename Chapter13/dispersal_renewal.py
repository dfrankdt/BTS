#!/usr/bin/env python3
"""
Dispersal Renewal: We construct the dispersion curve

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Nonlinearity
# =============================================================================
def f(u, A):
	y = u/(A + u)
	return y

# =============================================================================
# Main Simulation Function
# =============================================================================
def dispersal_renewal():
	# --- Parameters 
	bet = 2
	bsq = np.sqrt(bet**2 + 4)
	kp = (bsq + bet)/2
	km = (bsq - bet)/2

	# --- Variables
	L, Nx = 100, 2**9
	x = np.linspace(0, L, Nx + 1)
	XmY = np.outer( np.ones(Nx+1), x) - np.outer(x, np.ones(Nx+1))

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    dispersal_renewal()
