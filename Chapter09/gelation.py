#!/usr/bin/env python3
"""
Gelation: This script solves the pde

  dW/dt + d/dx (vW) = 0

where the initial profile is given by W0 and the velocity is prescribed
by a general form v = W/2 in order to give Burgers equation.

Produces: 

TO DO: Start with doGel
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# Nonlinearity (Initial Profile)
# =============================================================================
def W(x, C0, k):
	Wofz = C0*k*(x - x**(k-1))
	return Wofz

# =============================================================================
# Upwinding Velocity
# =============================================================================
def v(x, t, u):
	vel = 1/2*u
	return vel

# =============================================================================
# Upwinding Right-hand Side of the Differential Equation
# =============================================================================
def de_rhs(t, u, x):
	# --- u, x have Nx + 1 entries. 
	Nx = len(x) - 1
	dx = x[1] - x[0]
	
	# --- Need the flux J = v*u, but first need half step values for x and u 
	xjmh = np.zeros(Nx+2)
	xjmh[:-1] = x - dx/2
	xjmh[-1] = x[-1] + dx/2

	um = np.zeros(Nx+2)
	um[:-1] = u
	up = np.zeros(Nx+2)
	up[1:] = u
	
	ujmh = (um+up)/2

	# --- Call the velocity and get the flux, be sure to upwind.	
	vjmh = v(xjmh, t, ujmh)
	Jmh = vjmh*( (vjmh>0)*up + (vjmh<0)*um )
	
	# --- Now get the difference
	dJ = Jmh[1:] - Jmh[:-1]
	du = -dJ/dx
	
	return du
	
# =============================================================================
# Create Movie
# =============================================================================
def doMovie(x, t, U):
	# --- Initialize data structures
	Nt = len(t) - 1
	uinit = U[:,0]

	# --- Initialize movie
	fig, ax = plt.subplots()
	p_init = ax.plot(x, uinit, '--', label='Initial Profile')
	p_update = ax.plot([], [], '-r', label='Time Evolution')[0]
	ax.set(xlabel='z', ylabel='W(z, t)')
	ax.legend(loc='upper right')

	# --- Function to update the plot with the current frame
	def update(frame):
		tk = t[frame]
		uk = U[:, frame]
		p_update.set_xdata(x)
		p_update.set_ydata(uk)
		ax.set(title=f'Time t = {tk:.2f} s')
		return(p_update)

	ani = manimation.FuncAnimation(fig=fig, func=update, frames=range(Nt+1), interval=100)
	return ani
	
# =============================================================================
# Create Snapshots
# =============================================================================
def doSnapShots(x, t, U):
	# --- Initialize data structures
	Nt = len(t) - 1
	
	fig, ax = plt.subplots()
	u0 = U[:, 0]
	ax.plot(x, u0, '--')
	for kt in range(1, Nt+1):
		uk = U[:, kt]
		ax.plot(x, uk)
	ax.set(xlabel = 'z', ylabel = 'W(z, t)')
	ax.set(title = 'Numerical Solution via Upwinding')
	return fig, ax

# =============================================================================
# Exact Solution via Characteristics
# =============================================================================
def doExact(x, t, C0, k):
	# --- Initialize data structures
	Nt = len(t) - 1
	Nx = len(x) - 1
	w = np.zeros( (Nx+1, Nt+1) )
	z = np.zeros( Nx+1 )
	
	fig, ax = plt.subplots()
	w[:, 0] = W(x, C0, k)
	z[:] = x
	ax.plot(z, w[:, 0], '--')
	for kt in range(1, Nt+1):
		z[:] = W(x, C0, k)*t[kt] + x
		w[:, kt] = W(x, C0, k)
		ax.plot(z, w[:, kt])

	ax.set(xlabel = 'z', ylabel = 'W(z, t)')
	ax.set(title = 'Exact Solution via Characteristics')	
	return fig, ax	

# =============================================================================
# Characteristics via Resultant
# =============================================================================
def doChars(x, t, C0, k):
	# --- Initialize Plot
	fig, ax = plt.subplots()
	ax.set(xlabel = 'z', ylabel = 't')

	# --- First, characteristics via resultant constraints: eqns (9.95)
	w = W(x, C0, k)
	for kx in range(len(x)):
		ax.plot([x[kx], w[kx]*t[-1] + x[kx]], [0, t[-1]])
		

	# --- Use non-zero values of t, compute z via resultant
	tnz = t[1:]
	z = (9*C0**2*tnz**2 + 6*C0*tnz + 1)/(12*C0*tnz)
	ax.plot(z, tnz, '--', label = 'Envelope of Double Valued Solutions')
	ax.set(title = 'Characteristic Curves via Resultant')
	ax.legend(loc = 'upper right')
	return fig, ax

# =============================================================================
# Rootfinding
# =============================================================================
def getz0(C0, k, t):
	# --- Rootfinding to identify z0 for which z = 1
	xa = 0 
	fa = -1
	xb = 0.9999
	fb = W(xb, C0, k)*t + xb - 1
	
	for iter in range(20):
		xc = (xa+xb)/2
		fc = W(xc, C0, k)*t + xc - 1
		
		ftest = ((fc*fa) > 0)
		xa = ftest*xc + (1-ftest)*xa
		fa = ftest*fc + (1-ftest)*xa
		xb = (1-ftest)*xc + ftest*xb
		fb = (1-ftest)*fc + ftest*fb
	return xc

# =============================================================================
# DE RHS for Gelation
# =============================================================================
def de_rhsW(x, y, p):
	C0, k, tk = p
	z = W(x, C0, k)*tk + x
	w = W(x, C0, k)
	wp = C0*k*( 1 - (k-1)*x**(k-2))*tk + 1
	dw = w*wp
	return dw

# =============================================================================
# Gel Calculation: Monomer in gel vs monomer in polymer
# =============================================================================
def doGel(C0, k, tend):
	# --- Structures
	t = np.linspace(0, tend, 2**+1)
	W_Int = np.zeros(len(t))
	W1 = np.zeros(len(t))

	for kt in range(len(t)):
		xend = getz0(C0, k, t[kt])
		tspan = [0, xend]
		IVPW_pars = [C0, k, t[kt]]
		soln = solve_ivp(de_rhsW, tspan, [0], args = [IVPW_pars], t_eval=[xend])
		print(soln.y[0][0])

		W_Int[kt] = soln.y[0][0]
		W1[kt] = W(xend, C0, k)
	
	R = 3*C0/(3*C0*t + 1)
	fig, ax = plt.subplots()
	ax.plot(t, W1)
	
	return fig, ax
		
	
# =============================================================================
# Main Simulation Function
# =============================================================================
def gelation():
	# --- Discretizations
	L, Nx = 1, 100
	z = np.linspace(0, L, Nx+1)

	# --- Parameters
	C0 = 2/3
	k = 3

	# --- Initial Profile	
	w0 = W(z, C0, k)
	
	# --- Set the ODE for upwinding
	tf, Nt = 1.5, 15
	soln = solve_ivp(de_rhs, [0, tf], w0, args=[z], dense_output=True)

	# --- Structure to produce visualization
	t = np.linspace(0, tf, Nt+1)
	Wappx = soln.sol(t)
	#ani = doMovie(z, t, Wappx)
	fig, ax = doSnapShots(z, t, Wappx)
	
	# --- Exact Solution via Characteristics
	fig_ex, ax_ex = doExact(z, t, C0, k)
	
	# --- Characteristics via Resultant (just use a selection of the Nx+1 values)
	fig_char, ax_char = doChars(z[range(0, Nx+1, 10)], t, C0, k)
	
	# --- Gelation 
	tend = 2.5
	fig_gel, ax_gel = doGel(C0, k, tend)
	plt.show()

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    gelation()
