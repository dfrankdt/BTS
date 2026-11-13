#!/usr/bin/env python3
"""
Quorum Sensing Phase Planes:

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
	s0 = 0.05		# Basal production rate
	K = 1			# Production rate Kd
	g = 1.5			# Scale factor
	d = 1			# Transmembrane permeability
	alph = 0.0375	# Extracellular decay
	return s0, K, g, d, alph

# =============================================================================
# Nonlinearity
# =============================================================================
def f(u):
	s0, K, g, d, alph = doPars()
	y = s0 + g * u**2/(K**2 + u**2)
	return y

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(t, z, rho):
	s0, K, g, d, alph = doPars()
	u, v = z
	du = f(u) + d*(v - u)
	dv = -alph*v +  (rho*d)/(1 - rho)*(u - v) 
	dz = np.array([du, dv])
	return dz
		
# =============================================================================
# Phase Plane
# =============================================================================
def doPhasePlane(rho):
	# --- Initialization
	s0, K, g, d, alph = doPars()
	fig, ax = plt.subplots()
	ax.set(xlabel = 'u', ylabel = 'v', 	xlim = (0, 3), ylim = (0, 1.5))
	ax.set(title = rf'Phase Plane for $\rho$ = {rho:1.3f}')

	# --- Nullclines
	uu = np.linspace(0, 7, 2**9+1)
	ax.plot(uu, uu - f(uu)/d, '--k')
	ax.plot(uu, uu/((1-rho)/(rho) * alph/d + 1), '--k')
	
	# --- Solution trajectories
	tf = 100   # Large enough to achieve steady state
	v0_list = np.array([.15, .20, .25])
	for kv0 in range(len(v0_list)):
		v0 = v0_list[kv0]
		z0 = [0, v0]
		soln = solve_ivp(de_rhs, [0, tf], z0, args = [rho], dense_output=True)
		t = np.linspace(0, tf, 2**8+1)
		u, v = soln.sol(t)
		ax.plot(u, v)
	
	return fig, ax
		
# =============================================================================
# Main Simulation Function
# =============================================================================
def quorum_pp():
	# --- Phase planes for varying rho values
	rho = 0.007
	fig1, ax1 = doPhasePlane(rho)
	ax1.annotate('', xytext=(0.1, 0.3), xy=(0.3, 0.275), 
			arrowprops=dict(arrowstyle='->', lw=2))
	ax1.annotate('', xytext=(2, 0.5), xy=(1.7, 0.4),
			arrowprops=dict(arrowstyle='->', lw=2))
	
	rho = 0.02
	fig2, ax2 = doPhasePlane(rho)
	ax2.annotate('', xytext=(0.1, 0.3), xy=(0.3, 0.275),
			arrowprops=dict(arrowstyle='->', lw=2))
	ax2.annotate('', xytext = (.9, 0.25), xy = (1.25, 0.4),
			arrowprops=dict(arrowstyle='->', lw=2))

	rho = 0.035
	fig3, ax3 = doPhasePlane(rho)
	ax3.annotate('', xytext=(0.1, 0.3), xy=(0.3, 0.275),
			arrowprops=dict(arrowstyle='->', lw=2))
	ax3.annotate('', xytext=(1.5, 0.6), xy=(1.75, 0.75),
			arrowprops=dict(arrowstyle='->', lw=2))
			
	plt.show()		



# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    quorum_pp()
