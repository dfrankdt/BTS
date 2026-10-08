#!/usr/bin/env python3
"""
RBC_plots: Analysis to illustrate elements of a red blood cell production cycle.
The main result is a Hopf bifurcation, after which total number of cells in 
circulation becomes periodic. 

Produces:
 - Figure 1: Intersection of curves showing unique solution to equation (9.42), see
   Figure 9.1(a)
 - Figure 2: A decrease in cell death age gives a decrease in cell population
   coupled with an increase in the production of cells, see Figure 9.1(b)
 - Figure 3: Bifurcation diagram showing the Hopf bifurcation curve in the X/d
   (ratio of lifetime to delay) -- dA (product of rate and delay) plane, 
   see Figure 9.2(a)
 - Figure 4: Time dependent solution in region where solution is periodic,
   see Figure 9.2(b)
   
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# Nonlinearity
# =============================================================================
def F(N):
	# --- Hill Function
	y = 1/(1 + N**7)
	return y

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(t, z, p):
	# --- Identify parameters, state variables
	A, Nu, dx, NU, dy = p
	U = z[:NU+1]
	u = z[NU+1:]
	
	# --- Reprogram u0, U0
	u[0] = A*F(U[-1])
	U[0] = np.sum(u[0:Nu-1] + 2*u[1:Nu] + u[2:Nu+1])*(dx/2)
	
	# --- Identify fluxes
	Ju = -(u[1:Nu+1] - u[0:Nu])/dx
	JU = -(U[1:NU+1] - U[0:NU])/dy
	
	# --- Identify RHS
	dz = np.zeros( (NU+1) + (Nu+1) )
	dz[1:NU+1] = JU
	dz[NU+2:] = Ju
	return dz

# =============================================================================
# Figure 9.1
# =============================================================================
def do_Fig_9_1(b_values, N):
	# --- Figure (a)
	fig1, ax1 = plt.subplots()
	ax1.plot(N, F(N))
	ax1.set(xlabel = 'U', ylabel = r'Production Rate $f(U)/A$')
	for kb in range(len(b_values)):
		b = b_values[kb]
		ax1.plot(N, b*N, label = rf'$\beta = ${b:1.2f}')
	ax1.legend(loc='upper right')
	ax1.set(xlim = (-0.1, 2.1), ylim = (-0.1, 1.1))
	
	# --- Figure (b) 
	fig2, ax2 = plt.subplots()
	ax2.plot(N/F(N), N)
	ax2.set(xlabel = 'XA', ylabel = r'$U_0$')
	ax2.set(xlim = (-0.1, 4.1), ylim = (-0.1, 1.5))
	
	return fig1, fig2

# =============================================================================
# Figure 9.2
# =============================================================================
def do_Fig_9_2(d, X, A):
	# --- Figure (a)
	x = np.linspace(14/2**6, 14, 2**6+1)
	f = np.pi*x/( (2+x)*np.sin(2*np.pi/(2+x)) )
	p = 7
	N0p = f/(p-f)
	N0 = N0p**(1/p)
	dA = N0*(1+N0p)/x
	kp = np.where(N0p > 0)[0]
	
	fig1, ax1 = plt.subplots()
	ax1.set(xlabel = 'X/d', ylabel = 'dA')
	ax1.text(5, 0.15, 'Stable')
	ax1.text(5, 1.2, 'Unstable')
	ax1.plot(X/d, d*A, 'r.')
	ax1.plot(x[kp], dA[kp])
	ax1.set(xlim = (0, 14), ylim = (0, 2))
	
	# --- Figure (b)
	Nu = 2**7
	dx = X/Nu
	NU = 2**4
	dy = d/NU
	
	# --- IVP 
	tf = 500
	Uinit = np.ones(NU + 1)
	uinit = np.ones(Nu + 1)*dx/X
	
	zinit = np.zeros( (NU+1 + Nu+1) )
	zinit[:NU+1] = Uinit
	zinit[NU+1:] = uinit
	IVP_pars = [A, Nu, dx, NU, dy]
	
	soln = solve_ivp(de_rhs, [0, tf], zinit, args=[IVP_pars], dense_output=True)

	# --- Evaluate solution
	t = np.linspace(0, tf, 2**8+1)
	z = soln.sol(t)
	U = z[:NU+1, :]
	u = z[NU+1:, :]
	
	# --- Reprogram u0, U0
	u[0,:] = A*F(U[-1,:])
	U0 = np.sum(u[0:Nu-1,:] + 2*u[1:Nu,:] + u[2:Nu+1,:], 0)*(dx/2)
	fig2, ax2 = plt.subplots()
	ax2.plot(t, U0)
	ax2.set(xlabel = 'time (days)', ylabel = 'N(t)')
	ax2.set(ylim = (0,2.5))
	
	return fig1, fig2

# =============================================================================
# Main Simulation Function
# =============================================================================
def rbc_plots():
	# --- Steady State values for Figure 9.1
	b_values = np.array([0.8, 0.5, 0.2])
	N = np.linspace(0, 3, 2**8+1)

	# --- Do the plotting
	F91a, F91b = do_Fig_9_1(b_values, N)

	# --- Parameters for Figure 9.2
	d = 7		# Time Delay
	X = 50		# RBC lifetime
	A = 0.1		# production rate
	
	# --- Do the plotting
	F92a, F92b = do_Fig_9_2(d, X, A)
	
	plt.show()
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    rbc_plots()
