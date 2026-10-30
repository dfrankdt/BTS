#!/usr/bin/env python3
"""
Invasion: We construct the profile of invasion

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def invasion_profile():
	# --- Parameters
	beta = -2
	

	# --- Variables
	L, Nx = 5, 2**8
	x = np.linspace(-L, L, 2*Nx+1)


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    invasion_profile()
