#!/usr/bin/env python3
"""
Facilitated Diffusion

[Add some detail here]

Figures produced:
 - 
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Figure 1 (O2 Concentration)
# =============================================================================
def doFigOne(rho, s0, sL, sig):
	# --- Calculations
	mu = 1/( (1 + s0)*(1 + sL) )
	u =  sig/(1 + sig)
	sig0 = s0 + rho*s0/(1 + s0)
	fac = (1 + mu*rho)*(s0 - sL)
	xi = -(sig + rho*u - sig0)/fac
	stot = sig + u
	
	# --- Plotting
	fig1, ax1 = plt.subplots()
	ax1.plot(xi, sig, '-b', label = 'Free')
	ax1.plot(xi, u, '--r', label = 'Bound') 
	ax1.plot(xi, stot, '-.y', label = 'Total')
	ax1.set(xlabel = r'$\xi$', ylabel = 'Oxygen Concentration')
	ax1.set(xlim=(0, 1), ylim=(0, 3))
	ax1.legend(loc='upper right')

# =============================================================================
# Figure 2 (O2 Flux)
# =============================================================================
def doFigTwo(rho, s0, sL, sig):
	# --- Calculations
	mu = 1/( (1 + s0)*(1 + sL) )
	u =  sig/(1 + sig)
	sig0 = s0 + rho*s0/(1 + s0)
	fac = (1 + mu*rho)*(s0 - sL)
	xi = -(sig + rho*u - sig0)/fac
	sigp = fac*(1 + sig)**2/(sig**2 + rho + 2*sig + 1)
	up = sigp/(1 + sig)**2
	
	# --- Plotting
	fig2, ax2 = plt.subplots()
	ax2.plot(xi, sigp, '-b', label = 'Free Oxygen Flux')
	ax2.plot(xi, rho*up, '--r', label = 'Bound Oxygen Flux')
	ax2.set(xlabel = r'$\xi$', ylabel = 'Oxygen Flux')
	ax2.set(xlim = (0, 1), ylim = (0, 7))
	ax2.legend(loc = 'upper left')

# =============================================================================
# Figure 3 (External O2)
# =============================================================================
def doFigThree(rho, sig):
	# --- Calculations
	u =  sig/(1 + sig)
	gam = sig + rho*u

	# --- Plotting
	fig3, ax3 = plt.subplots()
	ax3.plot(gam, sig, '-b', label = r'$\rho = 5$')
	ax3.plot(sig, sig, '--r', label = r'$\rho = 0$') 
	ax3.set(xlabel = 'Oxygen Consumption', ylabel = 'External Oxygen Concentration')
	ax3.set(xlim = (0, 10), ylim = (0, 10))
	ax3.legend(loc = 'upper left')

# =============================================================================
# Figure 4 (Free O2)
# =============================================================================
def doFigFour(rho, sig, g):
	# --- Calculations
	u = sig/(1 + sig)
	gam = sig + rho*u
	
	xi = np.sqrt(gam/g)
	xi0 = np.sqrt(sig/g)

	s1 = g
	s0 = (-rho-s1-1+np.sqrt(4*rho*s1**2+rho**2+6*rho*s1+s1**2+2*rho+2*s1+1))/(1+s1)/2

	sig0 = np.linspace(s0, 20, 2**8+1)
	gam0 = sig0 + rho*sig0/(1 + sig0)
	xi1 = np.sqrt(abs(gam0 - s0 - rho*s0/(1 + s0))/g)
	
	# --- Plotting
	fig4, ax4 = plt.subplots()
	ax4.plot(xi, sig, '-b')
	ax4.plot(xi0, sig, '--r')
	ax4.plot(xi1, sig0, '-y')
	ax4.set(xlim=(0, 1), ylim=(0, 14))
	ax4.set(xlabel = 'Radius', ylabel = 'Free Oxygen')
	

# =============================================================================
# Main Simulation Function
# =============================================================================
def facilitated_diff_plots():
	# --- Global parameters
	rho = 15
	s0, sL = 2, 0.1
	sigma = np.linspace(0, 20, 2**8+1)

	# --- Create the figures
	doFigOne(rho, s0, sL, sigma)

	doFigTwo(rho, s0, sL, sigma)

	rho = 5
	doFigThree(rho, sigma)

	g = 14
	doFigFour(rho, sigma, g)

	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    facilitated_diff_plots()
