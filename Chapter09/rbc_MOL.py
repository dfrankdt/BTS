#!/usr/bin/env python3
"""
Method of Lines: We apply the method of lines to the PDE associated with RBC
production.

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
def F(z):
	y = 1/(1 + z**7)
	return y

# =============================================================================
# DE RHS
# =============================================================================
def de_rhs(t, z, p):
	A, Nu, dx, NU, dy = p
	U = z[:NU+1]
	u = z[NU+1:]
	
	u[0] = A*F(U[-1])
	U[0] = np.sum(u[0:Nu-1] + 2*u[1:Nu] + u[2:Nu+1])*(dx/2)
	
	Ju = -(u[1:Nu+1] - u[0:Nu])/dx
	JU = -(U[1:NU+1] - U[0:NU])/dy
	
	dz = np.zeros( (NU+1) + (Nu+1) )
	dz[1:NU+1] = JU
	dz[NU+2:] = Ju
	return dz

# =============================================================================
# Main Simulation Function
# =============================================================================
def rbc_MOL():
	# --- Discretizations
	d = 7
	X = 50
	A = 0.1
	tf = 500
	
	Nu = 2**7
	dx = X/Nu
	NU = 2**4
	dy = d/NU

	# --- IVP
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
	fig, ax = plt.subplots()
	ax.plot(t, U0)
	ax.set(xlabel = 'time (days)', ylabel = 'N(t)')
	ax.set(ylim = (0,2.5))
	
	plt.show()
	
	

# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
	rbc_MOL()
