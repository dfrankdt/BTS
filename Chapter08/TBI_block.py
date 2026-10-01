#!/usr/bin/env python3
"""
TBI Block: Identify the phase portrait and the critical relationship for the
blocked waves

Produces
 - Figure 1: Phase portrait in the U--W plane for the blocked wave, connects
   two integral curves (xi < 0 and xi > Y) to the solution of the DE in
   traveling wave coordinates, see Figure 8.8(a)
 - Figure 2: Critical relationship in the Y--A plane for existence of blocked
   waves, see Figure 8.8(b)

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
	# --- This event ensures w > 0
	u, w = z
	return w

w_zero.terminal = True
w_zero.direction = -1

def u_u1(x, z, p):
	# --- This event marks when the solution of the IVP has u = u1
	u0, u1, u2, A = p
	u, w = z
	return u-u1

u_u1.terminal = False
u_u1.direction = 1
	

# =============================================================================
# Block Trajectory
# =============================================================================
def getUY(zeros, Us, Ws):
	U0, U1, U2 = zeros

	Arstar = (F(U2,zeros) - F(U1,zeros) + F(Us,zeros))/(F(Us,zeros) - F(U1,zeros))
	IVP_args = [U0, U1, U2, Arstar]
	z0 = [Us, Ws]
	tmax = 100
	soln = solve_ivp(de_rhs, [0, tmax], z0, args = [IVP_args],
				events = [u_u1, w_zero], dense_output = True)
	Y = soln.t_events[0]
	tf = max(soln.t)
	t = np.linspace(0, tf, 2**15+1)
	z = soln.sol(t)
	UY, WY = soln.sol(Y)
	U, W = z
	return Y, U, W, UY, WY

# =============================================================================
# Critical Curve in Y - A
# =============================================================================
def getAYcrit(zeros):
	U0, U1, U2 = zeros
	
	# --- Discretize interval between U0 and U1 for initial points
	gamma = 0.9  # Don't get too close to U1
	NU = 2**8    # Number of points
	dU = (gamma*U1 - U0)/NU
	Us_values = np.linspace(U0+dU, gamma*U1-dU, NU-1)
	
	# --- Initialize Y and A
	Y = np.zeros(NU-1)
	A = np.zeros(NU-1)
	
	# --- Step through initial points
	for kU in range(len(Us_values)):
		Us = Us_values[kU]
		Ws = np.sqrt(-2*F(Us, zeros))

		Arstar = (F(U2,zeros) - F(U1,zeros) + F(Us,zeros))/(F(Us,zeros) - F(U1,zeros))
		IVP_args = [U0, U1, U2, Arstar]
		z0 = [Us, Ws]
		tmax = 100
		soln = solve_ivp(de_rhs, [0, tmax], z0, args = [IVP_args],
					events = [u_u1, w_zero], dense_output = True)
		Y[kU] = soln.t_events[0][0]
		A[kU] = Arstar
	
	return Y, A

# =============================================================================
# Fig 8.8 (a)
# =============================================================================
def getFig88a(zeros, Us, Ws):
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
	ax.plot(U1*np.ones(len(w)), w, '--')
	
	# --- Intermediate Trajectory satisfying Ar*U'' + f(U) = 0
	Y, U, W, UY, WY = getUY(zeros, Us, Ws)
	ax.plot(U, W, '--')
	ax.plot(UY, WY, 'ok')
	ax.plot(Us, Ws, 'ok')
	ax.text(0.05, 0, r'$(U_s, W_s)$')

	return fig, ax

# =============================================================================
# Fig 8.8 (b)
# =============================================================================
def getFig88b(alpha_vals, color_vals):
	fig, ax = plt.subplots()
	ax.set(xlim = (0, 50), ylim = (0, 50))
	ax.set(xlabel = 'Y', ylabel = r'$A_r$')

	for ka in range(len(alpha_vals)):
		alpha = alpha_vals[ka]
		U0, U1, U2 = 0, alpha, 1
		zeros = [U0, U1, U2]
		Y, A = getAYcrit(zeros)
		km = np.argmin(Y)
		ax.plot(Y[:km], A[:km], '-', color = color_vals[ka],
			label = rf'$\alpha$ = {alpha:1.2f}')
		ax.plot(Y[km:], A[km:], '--', color = color_vals[ka])
		ax.plot(Y[km], A[km], 'ok')

	ax.text(30, 35, 'blocked waves')
	ax.legend(loc = 'upper right')
	return fig, ax

# =============================================================================
# Main Simulation Function
# =============================================================================
def TBI_block():
	# --- Figure 8.8(a): Phase Portrait for alpha = 0.25 
	alpha = 0.25
	zeros = np.array([0, alpha , 1])
	Us = 0.05
	Ws = np.sqrt(-2*F(Us, zeros))
	fig88a, ax88a = getFig88a(zeros, Us, Ws)
	

	# --- Figure 8.8(b): Critical Y-A Curve
	alpha_list = [0.44, 0.3, 0.25]
	colorlist = np.array(['b', 'g', 'r'])
	fig88b, ax88b = getFig88b(alpha_list, colorlist)

	plt.show()	

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    TBI_block()
