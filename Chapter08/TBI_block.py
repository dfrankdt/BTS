#!/usr/bin/env python3
"""
TBI Block:

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
# Function Definitions
# =============================================================================
def f(x, zeros):
	# --- Cubic nonlinearity
	x0, x1, x2 = zeros
	y = (x-x0)*(x1-x)*(x-x2)
	return y

def F(x, a):
	# --- Integral of the cubic nonlinearity
	y = -(x**4/4 - (a+1)/3*x**3 + a/2*x**2)
	return y

# =============================================================================
# Main Simulation Function
# =============================================================================
def TBI_block():
	# --- Parameters
	Ar = 


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    TBI_block()
