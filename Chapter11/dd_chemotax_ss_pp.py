#!/usr/bin/env python3
"""
Steady State Analysis for Chemotaxis:

Note that we define kappa and compute the nontrivial equilibrium vstar.

Produces
 - Figure 1: Phase plane illustrating saddle, center, and periodic orbits, 
   see Figure 11.4(a)
 - Figure 2: Bifurcation diagram illustrating lengths for which solution exists,
   see Figure 11.4(b)
 - Figure 3: Typical solution trajectories, see Figure 11.5
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
# Nonlinearites
# =============================================================================
def f(v, k):
	# --- NL function f
	y = np.exp(k*(v - 1)) - v
	return y

def fp(v, k):
	# --- Derivative of f
	yp = k*np.exp(k*(v - 1)) - 1
	return yp

def fzero(k):
	# --- Newton's Method to ID the nontrivial zero of f
	zn = 0
	check = 1
	while check > 1e-14:
		znp1 = zn - f(zn, k)/fp(zn, k)
		check = np.abs(znp1 - zn)
		zn = znp1
	return zn

def F(v, k):
	# --- integral of f
	Y = 1/k*np.exp(k*(v - 1)) - v**2/2
	return Y	

def Fzero(v0, k):
	# --- Given (v0, 0) is on a trajectory, find (v1, 0)
	vn = 1.1
	check = 1
	
	while check > 1e-14:
		G = F(vn, k) - F(v0, k)
		Gp = f(vn, k)
		vnp1 = vn - G/Gp
		check = np.abs(vnp1 - vn)
		check = np.abs(F(vnp1, k) - F(v0, k))
		vn = vnp1
	return vn

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(xi, z, k):
	v, w = z
	dv = w
	dw = -f(v, k)
	dz = np.array([dv, dw])
	return dz

def w_zero(xi, z, k):
	v, w = z
	return w
w_zero.terminal = True
w_zero.direction = -1

# =============================================================================
# Phase Plane: Given kappa, return a phase plane
# =============================================================================
def doPP(k):
	# --- Identify the equilibrium and create a list of v-values
	vbar = fzero(k)
	
	# --- be sure v0list has one more element than colorlist
	v0list = np.linspace(vbar+0.01, 1, 6)
	colorlist = np.array(['b', 'g', 'r', 'c',  'm'])
	
	# --- Initialize plot and step through list
	fig, ax = plt.subplots()
	for kv in range(len(v0list)-1):
		# --- Initialize and solve IVP
		v0 = v0list[kv]
		soln = solve_ivp(de_rhs, [0, 10], [v0, 0], args = [k],
			 events = w_zero, dense_output = True)
		
		# --- Values for phase plane
		tf = max(soln.t)
		t = np.linspace(0, tf, 2**8+1)
		v, w = soln.sol(t)

		# --- Do the plotting
		ax.plot(v, w, color = colorlist[kv])
		ax.plot(v, -w, color = colorlist[kv])

		# --- Note the following should also plot the phase plane but 
		# --- perhaps there is some numerical error
		#v1 = Fzero(v0, k)
		#v = np.linspace(v0, v1, 2**3+1)
		#w = np.sqrt(2*(F(v0, k) - F(v, k)))
		#ax.plot(v, w, '.', color = colorlist[kv])
		#ax.plot(v, -w, '.', color = colorlist[kv])

	# --- Add equilibria to the phase plane
	ax.plot([vbar, 1], [0, 0], 'ok')
	ax.set(xlabel = 'v', ylabel = r'dv/d$\xi$')
	return fig, ax

# =============================================================================
# Bifurcation Diagram: Given kappa, return a bifurcation diagram
# =============================================================================
def doBif(k):
	# --- Identify the equilibrium and create a list of v-values
	vbar = fzero(k)
	v0list = np.linspace(vbar+0.01, 1 - 0.001, 2**6+1)

	vmin = np.zeros(len(v0list))
	vmax = np.zeros(len(v0list))
	L = np.zeros(len(v0list))

	for kv0 in range(len(v0list)):
		# --- Initialize and solve IVP
		v0 = v0list[kv0]
		soln = solve_ivp(de_rhs, [0, 10], [v0, 0], args = [k],
			events = w_zero, dense_output = True)
		# --- Values for Bifurcation Diagram
		L[kv0] = 2*max(soln.t)
		vmin[kv0] = v0
		vmax[kv0] = soln.y[0, -1]
	fig, ax = plt.subplots()
	ax.plot(L, vmax, label = r'$v_{\text{max}}$')
	ax.plot(L, vmin, label = r'$v_{\text{min}}$')
	
	ax.set(xlabel = 'Length l', ylabel = 'v')
	ax.set(xlim = (0, 16), ylim = (0, 1.6))
	ax.legend(loc = 'upper left')
	return fig, ax

# =============================================================================
# Typical Solution
# =============================================================================
def doSolution(k):
	vbar = fzero(k)
	v0 = 0.5
	
	# --- Solve IVP on half a periodic trajectory
	soln = solve_ivp(de_rhs, [0, 10], [v0, 0], args = [k],
			events = w_zero, dense_output = True)

	# --- Variables for plotting
	tf = max(soln.t)
	t = np.linspace(0, tf, 2**8+1)
	v = soln.sol(t)[0, :]
	
	# --- Plotting, leveraging symmetry
	fig, ax = plt.subplots()
	ax.plot(t, v, 'b', label = 'v')
	ax.plot(2*(tf-t/2), v, 'b')

	ax.plot(t, np.exp(k*(v-1)), 'r--', label = 'u/u*')
	ax.plot(2*(tf-t/2), np.exp(k*(v-1)), 'r--')

	ax.set(xlabel = r'$\xi$', ylabel = 'Solution', ylim = (0, 2))
	ax.legend(loc = 'upper right')

	return fig, ax

# =============================================================================
# Main Simulation Function
# =============================================================================
def dd_chemotax_ss_pp():
	# --- Parameters
	kappa = 2

	# --- Phase plane, bifurcation diagram, and solution
	fig1, ax1 = doPP(kappa)
	fig2, ax2 = doBif(kappa)
	fig3, ax3 = doSolution(kappa)
	
	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    dd_chemotax_ss_pp()
