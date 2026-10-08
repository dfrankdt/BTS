#!/usr/bin/env python3
"""
OU Process:

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation
rng = np.random.default_rng()

# =============================================================================
# Main Simulation Function
# =============================================================================
def OU_process():
	# --- Parameters
	L = 4
	d = 1
	K = 1
	
	# --- Discretizations
	Nx = 2**10
	Nt = 2**3
	



# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    OU_process()
