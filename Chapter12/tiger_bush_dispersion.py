#!/usr/bin/env python3
"""
Dispersion Curves for Tiger Bush Model:

"""

# =============================================================================
# Packages
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as manimation
rng = np.random.default_rng()

# =============================================================================
# Resultant
# =============================================================================
#def R(

# =============================================================================
# Dispersion Curve
# =============================================================================
def doDispersion(v, m, a_list):
	Na = len(a_list)
	Nw = 2**7
	# --- Variables
	wmax = 3
	w = np.linspace(0, wmax, Nw+1)
	lam = np.zeros( (2, Nw+1) , dtype=complex)

	# --- Initialize graph, step through a
	fig, ax = plt.subplots()
	ax.set(xlabel = r'$\omega$', ylabel = r'Re($\lambda$)')
	for ka in range(Na):
		a = a_list[ka]
		# --- Equilibria
		wbp = (a + np.sqrt(a**2 - 4*m**2))/2
		wbm = (a - np.sqrt(a**2 - 4*m**2))/2
		nbp = m/wbp
		nbm = m/wbm
		
		# --- Compute eigenvalues as a function of lambda
		for kw in range(Nw+1):
			fn = m - w[kw]**2
			fw = nbm**2
			gn = -2*m
			gw = 1 - nbm**2 + 1j*v*w[kw]
			J = np.array([ [fn, fw], [gn, gw] ])
			lam[:, kw] = np.linalg.eigvals(J)
		ax.plot(w, lam[0,:].real, label = f'a = {a:1d}')
		kwmax = np.argmax(lam[0,:].real)
		ax.plot(w[kwmax], lam[0, kwmax].real, '.k')
		print(f'for a = {a:1.1f}')
		print(f'we find lam = {lam[0, kwmax].real:.4f} {lam[0,kwmax].imag:+.4f}j')
		print(f'and wavelength = {2*np.pi/w[kwmax]:1.3f}')
	ax.set(ylim=(-1.5, 1.5))
	ax.plot([0, wmax], [0, 0], '--k')
	ax.legend(loc = 'upper right')
	return fig, ax
	
# =============================================================================
# Dispersion Curve
# =============================================================================
def doStability(v):
    # --- From trace of Jacobian
    m = np.linspace(2, 5, 2**8+1)
    a = m**2/np.sqrt(m-1)

	# --- Initialize stability diagram
    fig, ax = plt.subplots()
    ax.plot([0, 10], [0, 5])
    ax.plot(a, m, '--')
    
    # --- Get level curves
    gm = np.linspace(5/2**8, 10, 2**8)
    m = np.linspace(5/2**8, 5, 2**8)
    [GM, M] = np.meshgrid(gm, m)
    A = 2*M + GM
    Alpha = A/M
    N0 = (Alpha + np.sqrt(Alpha**2 - 4))/2
    N02 = N0**2
    B = 1 + N0**2
    
    # --- Resultant
    R = B**2*M**2*v**6 + 2*B**4*M*v**4 + 4*B**3*M**2*v**4
    R = R - 8*B**3*M**N02*v**4 - 4*B**2*M**3*v**4 + 8*B**2*M**2*N02*v**4 
    R = R + 4*B*M**3*N02*v**4 + B**6*v**2 + 8*B**5*M*v**2 - 8*B**4*M**2*v**2
    R = R - 40*B**4*M*N02*v**2 + 76*B**3*M**2*N02*v**2 - 40*B**2*M**3*N02*v**2 
    R = R - 32*B**2*M**2*N02**2*v**2 + 32*B*M**3*N02**2*v**2 + 4*M**4*N02**2*v**2
    R = R + 4*B**7 - 4*B**6*M - 24*B**5*M*N02
    R = R + 24*B**4*M**2*N02 + 48*B**3*M**2*N02**2 - 48*B**2*M**3*N02**2
    R = R - 32*B*M**3*N02**3 + 32*M**4*N02**3
    
    C = plt.contour(GM, M, R, [0])

    ax.text(3, 3, 'I')
    ax.text(6, 1, 'II')
    ax.text(9.5, 4.5, 'III')
#    ax.set(xlabel = 'a', ylabel = 'm', xlim = (0, 10), ylim = (0, 5))

    return fig, ax

# =============================================================================
# Main Simulation Function
# =============================================================================
def tiger_bush_dispersion():
	# --- Parameters
	v = 20
	m = 2
	a_list = [5, 6, 7, 8]
	
	fig1, ax1 = doDispersion(v, m, a_list)
	fig2, ax2 = doStability(v)
	plt.show()
	
		


# =============================================================================
# Execute the simulation if the script is run directly
# =============================================================================
if __name__ == "__main__":
    tiger_bush_dispersion()
