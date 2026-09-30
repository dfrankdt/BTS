#!/usr/bin/env python3
"""
TBI Block:

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# Function Definitions
# =============================================================================
def f(x, zeros):
	# --- Cubic nonlinearity
	x0, x1, x2 = zeros
	y = (x-x0)*(x1-x)*(x-x2)
	return y

def F(x, zeros):
	# --- Numerical integral of the cubic nonlinearity
	x0, x1, x2 = zeros
	dx = np.linspace(x0, x, 2**15+1)
	y = np.trapezoid(f(dx, zeros), dx)
	return y

def Fzero(zeros):
	# --- Zero of F
	x0, x1, x2 = zeros
	check = 1
	xn = (x1+x2)/2
	
	while check > 1e-8:
		xnp1 = xn - F(xn, zeros)/f(xn, zeros)
		check = np.abs(xnp1 - xn)
		xn = xnp1
	return xn

# =============================================================================
# IVP Structure
# =============================================================================
def de_rhs(x, z, p):
	u0, u1, u2, A = p
	zeros = np.array([u0, u1, u2])
	u, w = z
	du = w/A
	dw = -f(u, zeros)
	dz = np.array([du, dw])
	return dz

def w_zero(x, z, p):
	u, w = z
	return w

w_zero.terminal = True
w_zero.direction = -1


# =============================================================================
# Block Trajectory
# =============================================================================
def getUY(zeros, Us, Ws, A):
	U0, U1, U2 = zeros

	Arstar = (F(U2,zeros) - F(U1,zeros) - F(Us,zeros))/(F(Us,zeros) - F(U1,zeros))
	IVP_args = [U0, U1, U2, Arstar]
	z0 = [Us, Ws]
	tmax = 100
	soln = solve_ivp(de_rhs, [0, tmax], z0, args = [IVP_args],
				events = w_zero, dense_output = True)
	tf = max(abs(soln.t))
	t = np.linspace(0, tf, 2**15+1)
	z = soln.sol(t)
	U, W = z
	kmax = np.argmin( (U[:-1] - U1)*(U[1:] - U1) )
	print(kmax)
	kmax = np.argmax(W)
	print(kmax)
	return t, kmax, U, W

# =============================================================================
# Fig 8.8 (a)
# =============================================================================
def getFig88a(zeros, Us, Ws, A):
	U0, U1, U2 = zeros
	
	# --- Initiate Plot
	fig, ax = plt.subplots()
	ax.set(xlabel = 'U', ylabel = 'W')
	
	# --- Lower Trajectory satisfying U'' + f(U) = 0
	Umax = Fzero(zeros)
	u = np.linspace(0, Umax, 2**8+1)
	w = np.zeros(len(u))
	for ku in range(len(u)):
		w[ku] = np.sqrt( -2*F(u[ku], zeros) )
	ax.plot(u, w)
	
	# --- Upper Trajectory satisfying U'' + f(U) = 0
	u = np.linspace(0, U2, 2**8+1)
	w = np.zeros(len(u))
	for ku in range(len(u)):
		w[ku] = np.sqrt(2*(F(U2, zeros) - F(u[ku], zeros)))
	ax.plot(u, w)
	
	# --- Intermediate Trajectory satisfying Ar*U'' + f(U) = 0
	t, kmax, U, W = getUY(zeros, Us, Ws, A)
	ax.plot(U, W, '--')
	ax.plot(U[kmax], W[kmax], 'ok')	
		
	
	
	return fig, ax

# =============================================================================
# Main Simulation Function
# =============================================================================
def TBI_block():
	# --- Parameters
	Ar = 1
	alpha_list = [0.4, 0.3, 0.25]

	# --- Phase Portrait for alpha = 0.25
	alpha = 0.25
	zeros = np.array([0, alpha , 1])
	Umax = Fzero(zeros)

	Us = 0.05
	Ws = np.sqrt(-2*F(Us, zeros))
	fig88a, ax88a = getFig88a(zeros, Us, Ws, Ar)
	plt.show()

	# -- Step through alpha
	for ka in range(len(alpha_list)):
		alpha = alpha_list[ka]
		U0, U1, U2 = 0, alpha, 1
		zeros = [U0, U1, U2]


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    TBI_block()
