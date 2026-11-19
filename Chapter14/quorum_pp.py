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
# 
# =============================================================================
def gev(u, r):
	# --- Initialization
	s0, K, g, d, alph = doPars()

	cf = 1/((1 - r)*alph/(d*r) + 1) - 1
	y = cf*u + f(u)/d
	return y
		
# =============================================================================
# Triphasic Steady States
# =============================================================================
def getSScurves(u):
	# --- Initialization
	s0, K, g, d, alph = doPars()
	fig, ax = plt.subplots()
	ax.set(xlabel = r'Cell Density ($\rho$)', ylabel ='Intracellular Autoinducer (u)')
	ax.set(xlim = (0, 0.04), ylim = (0, 3))
	rho = alph*(d*u - f(u))/(d*f(u) - alph*f(u) + d*alph*u)
	ax.plot(rho, u)
	ax.text(0.020, 2.5, r'$c>0$')
	ax.text(0.005, 2.5, r'$c<0$')
	
	# --- Find the steady states in u
	rmin, rmax = np.min(rho), np.max(rho)
	rj = np.linspace(rmin, rmax, 2**8+1)
	rhok = np.zeros( len(rj) )

	# --- Step through rho values
	for jr in range(len(rj)):
		rk = 0
		r = rj[jr]
		geV = gev(u, r)
		Texp = np.zeros(len(rj))
		
#		# --- Find triphasic steady-states
#		sgc = geV[:-1]*geV[1:] 
#		idx = np.where(sgc<0)
#		if len(idx) == 3:
#			rk += 1
#			tst = (np.abs(geV[idx]) > np.abs(geV[idx+1]))
#			idx = idx + tst
#			Uss = u[idx]
#			cf = (1 - r)*alph/(d*r) + 1
#			Vss = Uss/cf
#			
#			# --- Identify the speed sign expressions
#			Ex1 = -r*d*(Vss - Uss)**2 / (2*(1 - r)) - alpha*Vss**2/2
#			DelU = (Uss[2] - Uss[0])/100
#			uss = np.linspace(Uss[0], Uss[2], 101)
#			Fus = f(uss)
#			Ex2 = r*(np.sum(Fus - Fus[0]/2 - Fus[2]/2)) * DelU*(1-r)
#			Texp[rk] = Ex1[2] -Ex1[1] + Ex2
#		
#			tcross = Texp[:-1]*Texp[1:]
#			tj = np.where(tcross<=0)
#	
#			tst = np.abs(Texp[tj]) > np.abs(Texp[tj+1])
#			rhocross = rhok[tj + tst]
	
#	ax.plot(rhocross*np.ones(len(u)), u)
	ax.plot(0.013*np.ones(len(u)), u, '--r')
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

	# --- Triphasic steady states in Intracellular Autoinducer
	u = np.linspace(0, 3, 2**8+1)
	fig4, ax4 = getSScurves(u)			
	plt.show()		



# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    quorum_pp()
