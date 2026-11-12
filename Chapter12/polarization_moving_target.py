#!/usr/bin/env python3
"""
Moving Target Polarization:

TO DO: Figure this out.
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
def polarization_moving_target():
	# --- Parameters
	L = 1
	Du = 0.0009
	Dv = 0.06
	del = 1
	K = 1
	gam = 0.5
	ampmx = 0.1
	tgrate = 0.05
	tgmx = 0.7
	
	# --- Variables
	u = np.linspace(0, 2, 2**8+1)
	f = gam * u**2/(1 + u**2)
	Vt = 4.2
	
	# --- 


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    polarization_moving_target()
