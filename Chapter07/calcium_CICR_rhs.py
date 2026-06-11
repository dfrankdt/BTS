#!/usr/bin/env python3
"""
CICR: Calcium induced calcium release

Bistability in CICR

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Main Simulation Function
# =============================================================================
def calcium_CICR_rhs():
	# --- Global Parameters
	Vp = 0.9
	kf = 1.1
	ct = 2
	K1 = 0.1
	K2 = 0.08
	gm = 5.5
	Kp = 0.1

	p_list = np.array([0.12, 0.15, 0.2])


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    calcium_CICR_rhs()
