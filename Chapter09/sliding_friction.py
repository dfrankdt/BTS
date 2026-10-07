#!/usr/bin/env python3
"""
Sliding Friction: We find dimensionless force as a function of dimensionless
velocity. Obtaining the dimensionless force requires solving three differential
equations and using their steady state values.

Produces for different stickiness parameters Kd
 - Figure 1(a): Force as a function of velocity, see Figure 9.12(a)
 - Figure 1(b): Velocity as a function of driving fluid velocity, see Figure 9.12(b)
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# DE Right Hand Side
# =============================================================================
def de_rhs(xi, z, V):
	n, N, f = z
	dn = - np.exp(xi)*n/V
	dN = n
	df = xi*n
	
	dz = np.array([dn, dN, df])
	return dz

# =============================================================================
# Main Simulation Function
# =============================================================================
def sliding_friction():
	# --- Parameters
	eta = 0.0005
	Kd_values = np.array([0.5, 1.0, 2.0])
	
	# --- Discretization: dimensionless velocity and desired output
	#     Note that values change quicly for small v, so we use a logspace
	Nv, vf = 2**6, 400
#	v = np.linspace(vf/Nv, vf, Nv)
	v = np.logspace(-1, np.log10(vf), Nv)
	Fc = np.zeros(Nv)
	
	# --- IVP Initialization
	s0 = [1, 0, 0]    # Initial condition of IVP
	tf = 10           # Large enough to ensure steady state
	
	# --- Initialize plots
	fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (12.8, 4.8))
	ax1.set(xlabel = 'Dimensionless Velocity', ylabel = 'Dimensionless Force')
	ax2.set(xlabel = 'Dimensionless Driving Fluid Velocity', ylabel = 'Dimensionless Velocity')
	
	# --- Step through values
	for kKd in range(len(Kd_values)):
		Kd = Kd_values[kKd]
		for kv in range(Nv):
			V = v[kv]
			
			# --- Solve IVP
			IVP_pars = V
			soln = solve_ivp(de_rhs, [0, tf], s0, args=[IVP_pars])
			
			# --- Get steady states
			z = soln.y
			N_inf = z[1, -1]
			f_inf = z[2, -1]
			
			# --- Pick off force and driving fluid velocity
			n0 = 1/(Kd*V + N_inf)
			Fc[kv] = n0*f_inf

		# --- Do the plotting
		vf = v + Fc/eta
		ax1.plot(v, Fc, label = rf'$K_d = ${Kd:1.1f}')
		ax2.plot(vf, v, label = rf'$K_d = ${Kd:1.1f}')

	ax1.legend(loc = 'upper right')
	ax2.legend(loc = 'upper right')		
	plt.show()
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    sliding_friction()
