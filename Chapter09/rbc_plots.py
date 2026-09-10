#!/usr/bin/env python3
"""
RBC_plots: We create plots that look like figure 9.1 (more soon). 

NOTE probably add some titles, fig 9.2 needs more work
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# Nonlinearity
# =============================================================================
def F(N):
	# --- Hill Function
	y = 1/(1 + N**7)
	return y

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
	print(kp)
	
	fig1, ax1 = plt.subplots()
	ax1.set(xlabel = 'X/d', ylabel = 'dA')
	ax1.annotate('Unstable', xytext = (5, 1.2), xy = (8, 1.2))
	ax1.annotate('Stable', xytext = (5, 0.2), xy = (8, 0.2))
	ax1.plot(X/d, d*A, 'r.')
	ax1.plot(x[kp], dA[kp])
	ax1.set(xlim = (0, 14), ylim = (0, 2))
	
	
	
	
	# --- Figure (b)
	Nn, dx = 2**7, X/Nn
	Nm = 2**4, d/Nm
	tspan = np.linspace(0, 500, 501)
	s0 = [np.ones(Nm, 1), np.ones(Nn, 1)*dx/X)
	soln = solve_ivp(de_rhs, tspan, s0)
	S = soln.sol(tspan)
	
	n0 = A/(1+s(Nm)**7)
	
	
	fig2, ax2 = plt.subplots()
	ax2.set(xlabel = 'time (days)', ylabel = 'N(t)')
	
	return fig1, fig2
	
	

# =============================================================================
# Main Simulation Function
# =============================================================================
def rbc_plots():
	# --- Parameters for Figure 9.1
	d = 70		# Time Delay
	X = 50		# RBC lifetime
	A = 1/50	# initial condition
	
	# --- Steady state values
	b_values = np.array([0.8, 0.5, 0.2])
	N = np.linspace(0, 3, 2**8+1)

	# --- Do the plotting
	F91a, F91b = do_Fig_9_1(b_values, N)

	# --- Parameters for Figure 9.2
	d = 7
	X = 50
	A = 0.1
	
	# --- Do the plotting
	F92a, F92b = do_Fig_9_2(d, X, A)
	
	plt.show()
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    rbc_plots()
