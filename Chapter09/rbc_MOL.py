#!/usr/bin/env python3
"""
Method of Lines: We apply the mehtod of lines to the PDE associated with RBC
production.

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
# Nonlinearity
# =============================================================================
def f(z):
	y = 1/(1 + z**7)
	return y

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(t, z, p):
	d, X, A, Nx, dx, Ny, dy = p
	

# =============================================================================
# Main Simulation Function
# =============================================================================
def rbc_MOL():
	# --- Discretizations
	d = 7
	X = 50
	A = 0.1
	
	Nx = 2**7
	dx = X/Nx
	Ny = 2**4
	dy = d/Ny


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	rbc_MOL()
