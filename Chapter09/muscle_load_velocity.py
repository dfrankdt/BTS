#!/usr/bin/env python3
"""
Muscle Load Velocity: We compute the steady state distributions for the 
fraction of crossbridges n(x) as a function of displacement for different
values of the velocity of the actin filament relative to the myosin filament.

Produces
 - Figure 1: Data and fit for original Hill model, see Figure 9.13
 - Figure 2: Steady state distributions, see Figure 9.16(a)
 - Figure 3: Update of Figure 1 with Huxley model, see Figure 9.16(b)

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Hill Fit
# =============================================================================
def doFit(g, v):
	"""
	The fit satisfies (g + a)v = b(p_0 - g), or av + bg + c = gv for unknown
	values a, b, and c = -b p0.  We construct a coefficient matrix A and RHS z.
	"""

	# --- Construct coefficient matrix and RHS for linear system
	A = np.zeros( (len(g), 3) )
	A[:,0] = v
	A[:,1] = g
	A[:,2] = np.ones(len(g))

	z = -g*v
	
	# --- Find the least squares solution (solve A'Ax = A'z)
	x = np.linalg.lstsq(A, z)[0]
	
	# --- Recover desired coeffs
	a = x[0]
	b = x[1]
	p0 = -x[2]/b

	# --- Plot the curve
	p = np.linspace(0, max(g), 2**8+1)
	vfit = b*(p0 - p)/(p+a)
	return p0, p, vfit
	
# =============================================================================
# Hill Plot
# =============================================================================
def doHillFig(g, v):
	# --- Initialize plot
	fig, ax = plt.subplots()
	ax.plot(g, v, 'o', label = 'Data')

	# --- Get fit, plot
	p0, p, vfit = doFit(g, v)
	ax.plot(p, vfit, label = 'Fit')
	ax.set(xlabel = 'Load, g', ylabel = 'Velocity of shortening, cm/s')
	ax.legend(loc = 'upper right')
	
	return fig, ax
# =============================================================================
# Main Simulation Function
# =============================================================================
def muscle_load_velocity():
	# --- Hill data 
	g = np.array([2.6,8.05,16.1,27,40.6,51.8,66])
	v = np.array([4.05,2.96,1.84,1.17,0.6,0.31,0])
	fig1, ax1 = doHillFig(g, v)

	# --- Huxley model discretizations and parameters
	xm = np.linspace(-3, 0, 2**7+1)
	xp = np.linspace(0, 0.999, 2**7+1)
	G2 = 3.919
	F1 = 13/16
	Vlist = np.array([0.0001, 0.6, 2, 4])
	colorlist = ['b', 'g', 'r', 'y']	
	labellist = [r'$V = 0$', r'$V = 0.15V_{max}$', r'$V = 0.5V_{max}$', r'$V = V_{max}$']

	# --- Initialize Plot for Fraction of Crossbridges
	fig2, ax2 = plt.subplots()
	ax2.plot([0, 0], [0, 1], '--k')
	ax2.set(xlabel = 'x/h', ylabel = 'n(x)')
	
	for jV in range(len(Vlist)):
		Vp = Vlist[jV]
		
		nI = F1*(1 - np.exp(-1/Vp))*np.exp(xm*G2/(2*Vp))
		nII = F1*(1 - np.exp((xp**2 - 1)/Vp))
		
		ax2.plot(xm, nI, '-', color = colorlist[jV], label = labellist[jV])
		ax2.plot(xp, nII, '-', color = colorlist[jV])
	ax2.legend(loc = 'upper left')
	
	# --- Initialize Plot for New load-velocity curve
	fig3, ax3 = plt.subplots()
	ax3.set(xlabel = 'Scaled Load', ylabel = 'Scaled Velocity')
	
	p0, p, vfit = doFit(g, v)
	ax3.plot(g/p0, v/5, 'ok', label = 'Data')
	ax3.plot(p/p0, vfit/5, '--k', label = 'Hill Fit')
	
	Vh = np.linspace(4/2**7, 4, 2**7)
	F = 1 - Vh*(1 - np.exp(-1/Vh))*(1 + Vh/(2*G2**2))
	ax3.plot(F, Vh/4, '-k', label = 'Huxley Fit')
	
	ax3.legend(loc = 'upper right')
	plt.show()
	
	

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	muscle_load_velocity()
