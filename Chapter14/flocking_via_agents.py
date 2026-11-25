#!/usr/bin/env python3
"""
Agent Based Flocking Behavior:

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
# Time Evolution
# =============================================================================
def doTrajectory(x0, s0, v, Nt, dt):
	# --- Initialization
	N = len(s0)
	X = np.zeros( (N, Nt+1) )
	X[:, 0] = x0

	# --- Update position and state
	Nxp, jxp = sum(s0), np.nonzero(s0)
	Nxm, jxm = N - Nxp, np.nonzero(1-s0)

	dx = s0*v*dt + (1-s0)*v*dt + np.sqrt(2*D*dt)*rng.normal(N)
		
	

# =============================================================================
# Main Simulation Function
# =============================================================================
def flocking_via_agents():
	# --- Parameters
	N = 10
	L = 1
	v = 1
	dt = 1e-2
	Nt = 1e3
	
	# --- Initial state, position
	s0 = np.round(rng.uniform(0, 1, N))
	x0 = np.linspace(2*np.pi/N, np.pi, N)

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    flocking_via_agents()
