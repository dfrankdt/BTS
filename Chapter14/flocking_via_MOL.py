#!/usr/bin/env python3
"""
Flocking Behavior: We use the method of lines to identify flocking behavior

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# Main Simulation Function
# =============================================================================
def flocking_via_MOL():
	# --- Parameters
	v = 0.5		# advective velocity
	D = 0.01	# diffusion coefficient
	Ut = 1.65	# population density (popn = Ut*L), threshold about Ut = 1

	# --- Discretization
	L, Nx = 20, 100
	dx = L/N
	x = np.linspace(0, L-dx, N)
	
	s0 = np.round(rng.uniform(0, 1, N))
	


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	flocking_via_MOL()
