#!/usr/bin/env python3
"""
Quorum Sensing: We examine the first order system governing the steady-state
solution on the bounded domain.

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# Parameters
# =============================================================================
def doPars():
	# --- Parameters 
	s0 = 0.05
	gamma = 1.5
	alph = 0.0375
	K = 1
	d = 1
	D = 1
	return s0, K, gamma, d, alph
# =============================================================================
# Main Simulation Function
# =============================================================================
def quorum_bounded_domain_ss_pp():



# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    quorum_bounded_domain_ss_pp()
