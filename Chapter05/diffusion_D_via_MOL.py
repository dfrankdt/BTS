#!/usr/bin/env python3
"""
Diffusion with Dirichlet BCs via the Method of Lines

We approximate the solution of the Diffusion Equation (5.37) -- (5.38) with
using either Dirichlet boundary conditions using the method of lines

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(t, u, alpha):
	Nu = np.size(u)
	sc = 2*np.ones(Nu)
	sc[0] = 1
	sc[-1] = 1
	sp = np.zeros(Nu)
	sp[1:] = u[:Nu-1]
	sm = np.zeros(Nu)
	sm[:-1] = u[1:Nu]
	d2u = -sc*u + sp + sm
	du = alpha*( d2u )
	return du
	

# =============================================================================
# Main Simulation Function
# =============================================================================
def diffusion_D_via_MOL():
	# --- Global parameters
	L = 1			# length of domain
	D = 1			# diffusion coefficient
	u0, uL = 0, 0	# boundary conditions
	Nx = 2**5		# number of partitions in x
	dt = 0.01
	tf = 0.5
	
	# --- Variables
	x = np.linspace(0, L, Nx + 1)
	dx = L/Nx
	alpha = D/dx**2


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    diffusion_D_via_MOL()
