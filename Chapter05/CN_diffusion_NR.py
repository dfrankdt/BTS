#!/usr/bin/env python3
"""
Diffusion with Neumann or Robin BCs via Crank-Nicolson

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_Diffusion_NR():
	# --- Parameters
	L = 1
	D = 1
	
	# --- Discretizations
	Nt, Nx = 2**4, 2**5
	dt, dx = 	
	# --- Boundary Conditions
	del0, U0 = 0, 0	# For Robin -Du_x + del0 u0 = U0
	delL, UL = 0, 0 # For Robin  Du_x + delL uN = UL

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_Diffusion_NR()
