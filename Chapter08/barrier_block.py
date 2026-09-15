#!/usr/bin/env python3
"""
Barrier Block: When f is a cubic nonlinearity we can patch the solution together
exactly. In the case that f is cubic-like, we need a numerical approach

Figures produced:
 - Figure 1: Trajectory (U-W plane), see Figure 8.6(a)
 - Figure 2: Trajectory (xi-U plane), see Figure 8.6(b)
 - Figure 3: Length of curve Y as a function of U(0), see Figure 8.7(a)
 - Figure 4: Critical blocking

TO DO: Definitely getting lost in the sauce on this one
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# =============================================================================
# Nonlinearities
# =============================================================================
def f(x, a):
	# --- Typical cubic nonlinearity
	y = x*(1 - x)*(x - a)
	return y

def F(x, a):
	# --- Integral of the cubic nonlinearity
	y = -(x**4/4 - (a+1)/3*x**3 + a/2*x**2)
	return y

def Fzero(a):
	# --- Zero of F
	xc = 2*(a+1)/3 - np.sqrt( 2*(2*a - 1)*(a - 2) )/3
	return xc

def getUY(a, k, Us, Ws):
	# --- Find intersection of trajectories via Newton's Method
	g0 = F(1, a) - 1/2*Ws**2 + k/2*Us**2

	# --- Initialize Newton's Method
	check = 1
	uk = Us
	
	# --- Iterate
	while check > 1e-9:
		g = k/2*uk**2 + F(uk, a)
		gp = 2*uk + f(uk, a)
		ukp1 = uk - (g-g0)/gp
		check = np.abs(ukp1 - uk)
		uk = ukp1
	wk = np.sqrt(2*(F(1, a) - F(uk, a)))
	return uk, wk
	
def getBlockTraj(a, k, Us, Ws):
	# --- Identify the trajectory satisfying U'' - kappa U = 0, U(0) = Us, W(0) = Ws
	UY, WY = getUY(a, k, Us, Ws)
	A0 = 1/2*(Us + Ws/np.sqrt(k))
	B0 = 1/2*(Us - Ws/np.sqrt(k))
	m = (UY + np.sqrt(UY**2 - 4*A0*B0))/(2*A0)
	Y = 1/np.sqrt(k)*np.log(m)
	
	xi = np.linspace(0, Y, 2**8+1)
	u = A0*np.exp(np.sqrt(k)*xi) + B0*np.exp(-np.sqrt(k)*xi)
	w = np.sqrt(k)*(A0*np.exp(np.sqrt(k)*xi) - B0*np.exp(-np.sqrt(k)*xi))
	return xi, u, w, Y

# =============================================================================
# IVP Structure
# =============================================================================
def de_rhs(x, z, a):
	u, w = z
	du = w
	dw = -f(u, a)
	dz = np.array([du, dw])
	return dz

def w_zero(x, z, a):
	u, w = z
	return w

w_zero.terminal = True
w_zero.direction = -1

	
# =============================================================================
# Figure 8.6 (a)
# =============================================================================
def getFig86a(alpha, kappa, Us, Ws):
	# --- Initiate Plot
	fig, ax = plt.subplots()
	ax.set(xlabel = 'U', ylabel = 'W')
	
	# --- Lower trajectory satisfying U'' + f(U) = 0
	Umax = Fzero(alpha)
	u = np.linspace(0, Umax, 2**8+1)
	w = np.sqrt(-2*F(u, alpha))
	ax.plot(u, w)
	
	# --- Upper trajectory satisfying U'' + f(U) = 0
	u = np.linspace(0, 1, 2**8+1)
	w = np.sqrt(2*(F(1, alpha) - F(u, alpha)))
	ax.plot(u, w)
	
	# --- Get intersection
	UY, WY = getUY(alpha, kappa, Us, Ws)
	ax.plot([Us, UY], [Ws, WY], 'ok')
	
	# --- Block trajectory
	xi, u, w, Y = getBlockTraj(alpha, kappa, Us, Ws)
	ax.plot(u, w, '--y')

	# --- Annotations
	ax.annotate('U(0)', xy = (Us-0.05, Ws-0.02))
	ax.annotate('U(Y)', xy = (UY+0.02, WY))
	return fig

# =============================================================================
# Figure 8.6 (b)
# =============================================================================
def getFig86b(alpha, kappa, Us, Ws):
	# --- Initiate IVP, plot
	IVP_args = [alpha]
	fig, ax = plt.subplots()
	ax.set(xlabel = r'$\xi$', ylabel = 'U')

	# --- Block trajectory
	xi, u, w, Y = getBlockTraj(alpha, kappa, Us, Ws)
	ax.plot(xi, u, '--y')
		
	# --- Block Parameters
	Uy = max(u)
	Wy = max(w)
	ax.plot([0, Y], [Us, Uy], 'ok')
	
	# --- Pre-block trajectory (xi < 0)
	tmax = 10
	zinit = [Us, Ws]
	tspan = [0, -tmax]
	soln = solve_ivp(de_rhs, tspan, zinit, args = IVP_args, 
			events = w_zero, dense_output = True)
	tf = max(abs(soln.t))
	t = np.linspace(0, -tf, 2**8+1)
	z = soln.sol(t)
	u = z[0, :]
	ax.plot(t, u)
	
	# --- Post-block trajectory (xi > Y)
	tmax = 20
	zinit = [Uy, Wy]
	tspan = [Y, Y+tmax]
	soln = solve_ivp(de_rhs, tspan, zinit, args = IVP_args, 
			events = w_zero, dense_output = True)
	tf = max(abs(soln.t))
	t = np.linspace(Y, tf, 2**8+1)
	z = soln.sol(t)
	u = z[0, :]
	ax.plot(t, u)
	
	# --- Annotations
	ax.annotate('(0, U(0))', xy = (1, Us))
	ax.annotate('(Y, U(Y))', xy = (Y+1, Uy))
	return fig

# =============================================================================
# Figure 8.7 (a)
# =============================================================================
def getFig87a(alpha, kappa, Us, Ws):
	# --- Initiate plot
	fig, ax = plt.subplots()
	ax.set(xlabel = 'U(0)', ylabel = 'Y')
	
	Umax = Fzero(alpha)
	Nu = 2**8
	du = Umax/Nu
	u0 = np.linspace(du, Umax, Nu)
	Y = np.zeros(len(u0))
	for ku in range(len(u0)):
		us = u0[ku]
		ws = np.sqrt(-2*F(us, alpha))
		xi, u, w, Y[ku] = getBlockTraj(alpha, kappa, us, ws)
	
	ax.plot(u0, Y)
	kmin = np.argmin(Y)
	ax.plot(u0[kmin], Y[kmin], '.k')
	
	ax.set(xlim = (0, 0.4), ylim = (0, 20))
	return fig
	
	
# =============================================================================
# Figure 8.7 (b)
# =============================================================================
def getFig87b(klist):
	# --- Initiate plot
	fig, ax = plt.subplots()
	ax.set(xlabel = r'$\alpha$', ylabel = 'Y')
	
	# --- Initiate Vars
	alist = np.linspace(0.5/2**8, 0.5, 2**8)
	Ycrit = np.zeros(len(alist))
	
	# --- Get critical lengths
	for kk in range(len(klist)):
		kappa = klist[kk]
		for ka in range(len(alist)):
			alpha = alist[ka]
			Umax = Fzero(alpha)
			u0 = np.linspace(Umax/2**6, Umax, 2**6)
			Y = np.zeros(len(u0))
			for ku in range(len(u0)):
				us = u0[ku]
				ws = np.sqrt(-2*F(us, alpha))
				xi, u, w, Y[ku] = getBlockTraj(alpha, kappa, us, ws)
			Ycrit[ka] = min(Y)
		ax.plot(alist, Ycrit, label = rf'$\kappa = ${kappa:1.3f}')
	ax.legend(loc = 'upper right')
	ax.set(ylim = (0, 20))
	ax.annotate('propagation failure', xy=(0.3, 10))
	ax.annotate('propagation success', xy = (0.05, 3))
	return fig

# =============================================================================
# Main Simulation Function
# =============================================================================
def barrier_block():
	# --- Parameters
	alpha = 0.245
	kappa = 0.05
	Us = 0.3
	Ws = np.sqrt(-2*F(Us, alpha))
	kappa_list = [0.1, 0.025, 0.01]

	fig86a = getFig86a(alpha, kappa, Us, Ws)
	fig86b = getFig86b(alpha, kappa, Us, Ws)	
	fig87a = getFig87a(alpha, kappa, Us, Ws)
	fig87b = getFig87b(kappa_list)
	plt.show()
	
# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    barrier_block()
