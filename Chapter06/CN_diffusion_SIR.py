#!/usr/bin/env python3
"""
CN Scheme to the RD Equation with SIR Dynamics

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

# =============================================================================
#  Nonlinearity
# =============================================================================
def F(u, v, alpha):
	# -- Reaction term for the RD Equation
	y = alpha * u * v
	return y

# =============================================================================
# Crank-Nicolson Method
# =============================================================================
def doCN():

# =============================================================================
# Create Movie
# =============================================================================
def doMovie():

# =============================================================================
# Main Simulation Function
# =============================================================================
def CN_diffusion_SIR():
	# --- Parameters
	D = 1
	r = 0.1
	eta = 0.5
	L = 80

	# --- Spatial and Temporal Scales
	tend = 1
	Nt, Nx = 2**4, 2**8
	dt, dx = tend/Nt,  L/Nx
	x = np.linspace(0, L, Nx+1)
	t = np.linspace(0, tend, Nt+1)

	# --- Initial Profile
	uinit = 0.1 * np.sech(2*(x - L/2))**2
	vinit = 1 - unit

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    CN_diffusion_SIR()
