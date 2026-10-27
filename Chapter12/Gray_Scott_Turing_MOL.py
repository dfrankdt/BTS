#!/usr/bin/env python3
"""
Gray-Scott-Turing: We examine neutral stability of the Turing Mechanism, 
the phase portrait for the Gray Scott model, and steady-state solutions 
of the Gray-Scott model via the Method of Lines.

TO DO: The MOL does not pick up the instability based on an initial condition
that is a perturbation of the steady state.
"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation

# =============================================================================
# DE Structure (ODE and PDE)
# =============================================================================
def ode_rhs(t, z, p):
	kv, rho = p
	u, v = z
	du = -u*v**2 + rho*(1 - u)
	dv = u*v**2 - (rho + kv)*v
	dz = np.array([du, dv])
	return dz

def pde_rhs(t, z, p):
	# --- Parameters and structure
	kv, rho, Du, Dv, dx = p
	N = int( (len(z) - 2)/2 )
	U = z[:N+1]
	V = z[N+1:]
	
	# --- Second difference operators with no flux at boundaries
	D2U = np.zeros(N+1)
	D2U[0] = -2*U[0] + 2*U[1]
	D2U[-1] = 2*U[-2] - 2*U[-1]
	D2U[1:N] = U[0:N-1] - 2*U[1:N] + U[2:N+1]
	D2U = D2U/dx**2
	
	D2V = np.zeros(N+1)
	D2V[0] = -2*V[0] + 2*V[1]
	D2V[-1] = 2*V[-2] - 2*V[-1]
	D2V[1:N] = V[0:N-1] - 2*V[1:N] + V[2:N+1]
	D2V = D2V/dx**2
	
	# --- Right-hand side of the PDE
	F = U*V**2
	dUdt = Du*D2U - F + rho*(1 - U)
	dVdt = Dv*D2V + F - (rho + kv)*V
	
	# --- Evaluate for dz/dt
	dz = np.zeros( (N+1) + (N+1) )
	dz[:N+1] = dUdt
	dz[N+1:] = dVdt
	return dz
		
# =============================================================================
# Parameters
# =============================================================================
def doPars(kv, rho):
	# --- Computed Values
	p = (rho + kv)**2/rho
	u0 = (1 - np.sqrt(1 - 4*p))/2
	v0 = np.sqrt(rho*(1 - u0)/u0)

	fu = -v0**2 - rho
	fv = -2*u0*v0
	gu = v0**2
	gv= 2*u0*v0 - rho - kv
	det = fu*gv - fv*gu
	
	return p, u0, v0, fu, fv, gu, gv, det

# =============================================================================
# Neutral Stability
# =============================================================================
def doStability(kv, rho, Du, Dv, L):
	# --- Initialize Plot
	fig, ax = plt.subplots()
	ax.set(xlabel = r'$\omega^2 D_v$', ylabel = r'$\omega^2 D_u$')

	# --- Modes
	n = np.linspace(1, 16, 16)
	wn = n*np.pi/L
	ax.plot(wn**2*Dv, wn**2*Du, '.k')
	
	# --- Initialize Variables, parameters for stability curve
	p, u0, v0, fu, fv, gu, gv, det = doPars(kv, rho)
	dv = np.linspace(0, 0.08, 2**7+1)
	du = (det - dv*fu)/(gv - dv)
	
	# --- Values for stability curve
	a = gv**2
	b = 2*fu*gv - 4*det
	c = fu**2
	alpha = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
	
	# --- Graph
	ax.plot(dv, du, label = 'Neutral Stability')
	ax.plot(dv, alpha*dv, '--', label = r'$D_u = \alpha D_v$')
	ax.set(xlim = (0, gv), ylim = (0, 4))
	ax.text(0.02, 2, 'Unstable')
	ax.text(0.06, 1, 'Stable')
	ax.legend(loc = 'upper left')
	
	return fig, ax
# =============================================================================
# Phase Portrait
# =============================================================================
def doPP(kv, rho):
	# --- Initialize Plot
	fig, ax = plt.subplots()

	# --- Initialize Variables
	Nu = 2**10
	u = np.linspace(1/Nu, 1, Nu)
	p, u0, v0, fu, fv, gu, gv, det = doPars(kv, rho)	

	# --- Nullclines, EP for plotting
	v1 = np.sqrt(rho*(1 - u)/u)
	v2 = (rho + kv)/u
	us1 = (1 + np.sqrt(1 - 4*p))/2
	vs1 = np.sqrt(rho*(1 - us1)/us1)
	
	# --- ODE Solution
	z0 = [1, 0.095]
	tf = 50
	ODE_pars = [kv, rho]
	soln = solve_ivp(ode_rhs, [0, tf], z0, args = [ODE_pars], dense_output = True)
	t = np.linspace(0, tf, 2**8+1)
	U, V = soln.sol(t)
	
	# --- Do the plotting
	ax.plot(u, v1, '--', label = 'u nullcline')
	ax.plot(u, v2, '--', label = 'v nullcline')
	ax.plot([u0, us1], [v0, vs1], '.k')
	ax.plot(U, V, label = 'Trajectory')
	
	ax.set(xlim = (-0.1 , 1.1), ylim = (0, 1))
	ax.annotate("", xytext=(0.7, 0.4), xy=(0.5, 0.6), 
		arrowprops = dict(arrowstyle="->"))
	ax.legend(loc = 'upper right')
	return fig, ax

# =============================================================================
# Method of Lines
# =============================================================================
def doMOL(rho, kv, Du, Dv, L, n0):
	# --- Parameters
	p, u0, v0, fu, fv, gu, gv, det = doPars(kv, rho)

	# --- Spatial structure
	Nx = 2**8
	x = np.linspace(0, L, Nx+1)
	U0 = u0 - 0.02*np.cos(n0*np.pi*x/L)
	V0 = v0 - 0.02*np.cos(n0*np.pi*x/L)

	# --- MOL Structure
	z0 = np.zeros( (Nx + 1) + (Nx + 1))
	z0[:Nx+1] = U0
	z0[Nx+1:] = V0
	PDE_pars = [kv, rho, Du, Dv, Nx/L]
	tf = 2500

	# --- Use MOL to solve the PDE, we want the steady state
	soln = solve_ivp(pde_rhs, [0, tf], z0, args = [PDE_pars], dense_output = True)
	Uss = soln.y[:Nx+1,-1]
	Vss = soln.y[Nx+1:,-1]
	
	fig, (ax1,ax2) = plt.subplots(1, 2, figsize = (12.8, 4.8))
	ax1.plot(x, Uss, label = 'Steady State')
	ax1.plot(x, U0, '--', label = 'Initial Condition')
	ax1.set(xlabel = 'x', ylabel = 'u(x, t)', ylim = (0, 1))
	ax1.legend(loc = 'upper right')
	
	ax2.plot(x, Vss, label = 'Steady State')
	ax2.plot(x, V0, '--', label = 'Initial Condition')
	ax2.set(xlabel = 'x', ylabel = 'v(x, t)', ylim = (0, 1))
	ax2.legend(loc = 'upper right')

	return fig, (ax1, ax2)
# =============================================================================
# Main Simulation Function
# =============================================================================
def Gray_Scott_Turing_MOL():
	# --- Parameters
	rho = 1/16
	kv = 0.02
	Du = 1
	Dv = 0.025
	L = 20
		
	# --- Neutral Stability Curve and Phase Plane
	fig1, ax1 = doStability(kv, rho, Du, Dv, L)
	fig2, ax2 = doPP(kv, rho)

	# --- Structure for Method of Lines
	n0 = 6
	fig3, (ax3a, ax3b) = doMOL(rho, kv, Du, Dv, L, n0)

	n0 = 8
	fig4, (ax4a, ax4b) = doMOL(rho, kv, Du, Dv, L, n0)

	plt.show()


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    Gray_Scott_Turing_MOL()
