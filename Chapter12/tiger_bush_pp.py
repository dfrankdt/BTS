#!/usr/bin/env python3
"""
Phase Plane for the Tiger Bush System:

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
# DE Structure
# =============================================================================
def de_rhs(t, z, p):
	a, m = p
	n, w = z
	dn = w*n**2 - m*n
	dw = a - w - w*n**2
	dz = np.array([dn, dw])
	return dz

# =============================================================================
# Phase Plane
# =============================================================================
def getPP(a, m):
	# --- Initialize Figure
	fig, ax = plt.subplots()
	ax.set(xlabel = 'n', ylabel = 'w', xlim=(-0.1, 4.1), ylim = (-0.1, 6.1))
	
	# --- Nullclines
	x = np.linspace(4/2**8, 4, 2**8)
	yw = a/(1+x**2)
	yn = m/x
	ax.plot(x, yn, '--')
	ax.plot(x, yw, '--')
	
	# --- Solution Trajectories
	tf = 10
	t = np.linspace(0, tf, 2**8+1)
	IVP_pars = [a, m]

	z0 = [1, 0.]
	soln = solve_ivp(de_rhs, [0, tf], z0, args = [IVP_pars], dense_output=True)
	z = soln.sol(t)
	n, w = z
	ax.plot(n, w, 'k')
	
	z0 = [1.5, 0.]
	soln = solve_ivp(de_rhs, [0, tf], z0, args = [IVP_pars], dense_output=True)
	z = soln.sol(t)
	n, w = z
	ax.plot(n, w, 'k')

	# --- Equilibrium Solutions: -wb**2 + a*wb - m**2 = 0, nb = m/wb
	wbp = (a + np.sqrt(a**2 - 4*m**2))/2
	wbm = (a - np.sqrt(a**2 - 4*m**2))/2
	nbp = m/wbp
	nbm = m/wbm
	ax.plot([nbm, nbp], [wbm, wbp], '.k')

	return fig, ax

# =============================================================================
# Stability
# =============================================================================
def getBifDiag():
	# --- From trace of Jacobian
	m = np.linspace(2, 5, 2**8+1)
	a = m**2/np.sqrt(m-1)
	
	fig, ax = plt.subplots()
	ax.plot([0, 10], [0, 5])
	ax.plot(a, m, '--')
	ax.text(3, 3, 'I')
	ax.text(6, 1, 'II')
	ax.text(9.5, 4.5, 'III')
	ax.set(xlabel = 'a', ylabel = 'm', xlim = (0, 10), ylim = (0, 5))
	
	return fig, ax

# =============================================================================
# Main Simulation Function
# =============================================================================
def tiger_bush_pp():
	# --- Parameters
	a = 6
	m = 2
	
	# --- Phase Plane
	fig1, ax1 = getPP(a, m)
	fig2, ax2 = getBifDiag()

	plt.show()
	
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	tiger_bush_pp()
